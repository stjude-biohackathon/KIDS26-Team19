"""Shared schema and streaming header parser for the GEO DuckDB.

Used by:
  parsing/initDb.py    creates an empty database with the fixed tables
  parsing/updateDb.py  adds series matrix files to an existing database

Tables: diagnosis -> dataset (one row per series matrix file) -> sample (one row
per GEO sample). Treatment and control are derived from all !Sample_characteristics_ch1
header lines (treatment/drug row), not stored as raw characteristics.
"""

from __future__ import annotations

import csv
import gzip
import io
import os
import re

TABLE_BEGIN_MARKER = "!series_matrix_table_begin"

DEFAULT_DB_PATH = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "geo.db")
)

WANTED_KEYS = (
    "Series_geo_accession",
    "Series_platform_id",
    "Sample_geo_accession",
    "Sample_organism_ch1",
    "Sample_data_row_count",
    "Sample_molecule_ch1",
    "Series_pubmed_id",
)

HEADER_SAMPLE_KEYS = (
    "Sample_geo_accession",
    "Sample_organism_ch1",
    "Sample_data_row_count",
    "Sample_molecule_ch1",
)
DERIVED_SAMPLE_COLUMNS = ("treatment", "control", "drug")

SERIES_KEYS = tuple(k for k in WANTED_KEYS if k.startswith("Series_"))

OPTIONAL_KEYS = ("Series_pubmed_id",)

KEY_COLUMNS = {key: key.lower() for key in WANTED_KEYS}
for col in DERIVED_SAMPLE_COLUMNS:
    KEY_COLUMNS[col] = col

INTEGER_COLUMNS = frozenset({"sample_data_row_count"})

DATASET_COLUMNS = ("diagnosis_id", "source_file", *(KEY_COLUMNS[k] for k in SERIES_KEYS))
SAMPLE_COLUMNS = (
    "dataset_id",
    "diagnosis_id",
    *(KEY_COLUMNS[k] for k in HEADER_SAMPLE_KEYS),
    *DERIVED_SAMPLE_COLUMNS,
)


def _column_type(column: str) -> str:
    return "BIGINT" if column in INTEGER_COLUMNS else "VARCHAR"


_sample_ddl_parts = [
    f"{KEY_COLUMNS[key]} {_column_type(KEY_COLUMNS[key])}" for key in HEADER_SAMPLE_KEYS
] + [f"{col} VARCHAR" for col in DERIVED_SAMPLE_COLUMNS]

_DDL = f"""
CREATE SEQUENCE IF NOT EXISTS diagnosis_id_seq START 1;
CREATE SEQUENCE IF NOT EXISTS dataset_id_seq START 1;
CREATE SEQUENCE IF NOT EXISTS sample_id_seq START 1;

CREATE TABLE IF NOT EXISTS diagnosis (
    diagnosis_id BIGINT PRIMARY KEY DEFAULT nextval('diagnosis_id_seq'),
    diagnosis_name VARCHAR NOT NULL UNIQUE
);

CREATE TABLE IF NOT EXISTS dataset (
    dataset_id BIGINT PRIMARY KEY DEFAULT nextval('dataset_id_seq'),
    diagnosis_id BIGINT NOT NULL REFERENCES diagnosis(diagnosis_id),
    source_file VARCHAR NOT NULL UNIQUE,
    {", ".join(
        f"{KEY_COLUMNS[key]} {_column_type(KEY_COLUMNS[key])}" for key in SERIES_KEYS
    )}
);

CREATE TABLE IF NOT EXISTS sample (
    sample_id BIGINT PRIMARY KEY DEFAULT nextval('sample_id_seq'),
    dataset_id BIGINT NOT NULL REFERENCES dataset(dataset_id),
    diagnosis_id BIGINT NOT NULL REFERENCES diagnosis(diagnosis_id),
    {", ".join(_sample_ddl_parts)}
);
"""

INSERT_DATASET_SQL = (
    "INSERT INTO dataset (dataset_id, {columns}) "
    "VALUES (nextval('dataset_id_seq'), {placeholders}) RETURNING dataset_id"
).format(
    columns=", ".join(DATASET_COLUMNS),
    placeholders=", ".join("?" for _ in DATASET_COLUMNS),
)

INSERT_SAMPLE_SQL = (
    "INSERT INTO sample (sample_id, {columns}) "
    "VALUES (nextval('sample_id_seq'), {placeholders})"
).format(
    columns=", ".join(SAMPLE_COLUMNS),
    placeholders=", ".join("?" for _ in SAMPLE_COLUMNS),
)


class MatrixError(Exception):
    """A series matrix file cannot be loaded."""


# --- Treatment / control from Sample_characteristics_ch1 ---

_EXCLUDED_TREATMENT_VALUES = frozenset({"untreated", "vehicle", "dmso", "control"})


def _strip_geo_quotes(text):
    if text is None:
        return ""
    s = str(text).strip()
    if len(s) >= 2 and s[0] == '"' and s[-1] == '"':
        s = s[1:-1]
    return s.strip()


def _normalize_for_match(text):
    s = _strip_geo_quotes(text).lower()
    s = re.sub(r"[^a-z0-9\s]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def _normalize_treatment_value(text):
    return _normalize_for_match(text)


def _is_excluded_treatment_value(value):
    if not value:
        return False
    norm = _normalize_treatment_value(value)
    for token in _EXCLUDED_TREATMENT_VALUES:
        if norm == token:
            return True
        if norm.startswith(token + " "):
            return True
        if token == "control" and norm.endswith(" control"):
            return True
    return False


def _cell_keyword(text):
    norm = _normalize_for_match(text)
    if "treatment" in norm:
        return "treatment"
    if "drug" in norm:
        return "drug"
    return None


def _label_mentions_treatment_or_drug(label_text):
    norm = _normalize_for_match(label_text)
    return "treatment" in norm or "drug" in norm


def _parse_labeled_treatment_cell(text):
    raw = _strip_geo_quotes(text)
    if not raw or not _cell_keyword(raw):
        return None
    value = None
    for sep in (":", "-", "/"):
        if sep in raw:
            left, right = raw.split(sep, 1)
            if _label_mentions_treatment_or_drug(left):
                value = _strip_geo_quotes(right)
                break
    if value is None:
        value = raw
    if not value or not str(value).strip():
        return None
    return value


def _is_control_like_value(value):
    norm = _normalize_treatment_value(value)
    if norm == "none":
        return True
    if "dmso" in norm:
        return True
    if "untreat" in norm:
        return True
    if "vehicle" in norm:
        return True
    if "control" in norm:
        return True
    if "ctrl" in norm:
        return True
    return _is_excluded_treatment_value(value)


def _extract_value_from_cell(text):
    value = _parse_labeled_treatment_cell(text)
    if value is None or _is_control_like_value(value):
        return None
    return value


def _extract_control_from_cell(text):
    value = _parse_labeled_treatment_cell(text)
    if value is None or not _is_control_like_value(value):
        return None
    return value


def _row_keyword(row):
    for cell in row:
        if _cell_keyword(cell) == "treatment":
            return "treatment"
    for cell in row:
        if _cell_keyword(cell) == "drug":
            return "drug"
    return None


def _select_characteristics_row(characteristics_rows):
    treatment_row = None
    drug_row = None
    for row in characteristics_rows:
        kw = _row_keyword(row)
        if kw == "treatment" and treatment_row is None:
            treatment_row = row
        elif kw == "drug" and drug_row is None:
            drug_row = row
    return treatment_row or drug_row


def extract_treatment_from_characteristics(characteristics_rows, n_samples):
    """Pick a Sample_characteristics_ch1 row mentioning treatment/drug; parse per sample."""
    selected = _select_characteristics_row(characteristics_rows)
    if selected is None:
        return [None] * n_samples
    out = []
    for i in range(n_samples):
        cell = selected[i] if i < len(selected) else None
        if cell is None:
            out.append(None)
            continue
        value = _extract_value_from_cell(cell)
        out.append(value if value else None)
    return out


def extract_control_from_characteristics(characteristics_rows, n_samples):
    """DMSO, vehicle, untreated, none, and control-like values from the treatment/drug row."""
    selected = _select_characteristics_row(characteristics_rows)
    if selected is None:
        return [None] * n_samples
    out = []
    for i in range(n_samples):
        cell = selected[i] if i < len(selected) else None
        if cell is None:
            out.append(None)
            continue
        value = _extract_control_from_cell(cell)
        out.append(value if value else None)
    return out


def derive_treatment_control_from_cell(cell_value: str | None) -> tuple[str | None, str | None]:
    """Best-effort derive treatment/control from one stored characteristics cell."""
    if cell_value is None or not str(cell_value).strip():
        return None, None
    rows = [[cell_value]]
    return (
        extract_treatment_from_characteristics(rows, 1)[0],
        extract_control_from_characteristics(rows, 1)[0],
    )


_DOSAGE_UNIT = r"\d+(?:\.\d+)?\s*(?:µM|μM|uM|nM|pM|mM)\s+"
_DOSAGE_PREFIX = re.compile(
    r"^\s*" + _DOSAGE_UNIT + r"(.+)$",
    re.IGNORECASE,
)
_DOSAGE_INLINE = re.compile(_DOSAGE_UNIT, re.IGNORECASE)
_DRUG_TOKEN_AFTER_DOSAGE = re.compile(r"([^,+]+)")
_DOSAGE_LEADING = re.compile(r"^\s*" + _DOSAGE_UNIT, re.IGNORECASE)
_INHIBITOR_IN_TEXT = re.compile(r"(?i)(.+?)\s+inhibitor\b")
_DURATION_PREFIX = re.compile(
    r"^\s*\d+(?:\.\d+)?\s*(?:hr|hrs|h|min|mins|minutes|day|days)\s+",
    re.IGNORECASE,
)


def _strip_leading_duration(text: str) -> str:
    text = text.strip()
    while text:
        match = _DURATION_PREFIX.match(text)
        if not match:
            break
        text = text[match.end() :].strip()
    return text


def _is_duration_only_prefix(prefix: str) -> bool:
    return not _strip_leading_duration(prefix.strip())


def _strip_leading_dosage(name: str) -> str:
    text = name.strip()
    while text:
        match = _DOSAGE_LEADING.match(text)
        if not match:
            break
        text = text[match.end() :].strip()
    return text


def _inhibitor_drug_labels(text: str) -> list[str]:
    stripped = text.strip()
    if not stripped:
        return []
    labels: list[str] = []
    for match in _INHIBITOR_IN_TEXT.finditer(stripped):
        target = match.group(1).strip().lstrip("+").strip()
        target = _strip_leading_dosage(target)
        if not target:
            continue
        label = f"{target} inhibitor"
        if label not in labels:
            labels.append(label)
    return labels


def _drug_from_inhibitor_suffix(text: str) -> str | None:
    labels = _inhibitor_drug_labels(text)
    return labels[0] if labels else None


def _merge_drug_parts(parts: list[str], extra: str | None) -> list[str]:
    if not extra or extra in parts:
        return parts
    return parts + [extra]


def _label_part_drug(part: str) -> str:
    part = part.strip()
    if not part:
        return part
    labels = _inhibitor_drug_labels(part)
    if labels:
        return labels[0]
    return part


def _collect_drugs_from_segment(segment: str) -> list[str]:
    """Split a post-leading-dosage name segment into one or more drug names."""
    segment = segment.strip()
    if not segment:
        return []
    match = _DOSAGE_INLINE.search(segment)
    if not match:
        segment = segment.rstrip("+").strip()
        if "+" in segment:
            parts = [part.strip() for part in re.split(r"\s*\+\s*", segment)]
            return [_label_part_drug(part) for part in parts if part.strip()]
        labeled = _label_part_drug(segment)
        return [labeled] if labeled else []
    prefix = segment[: match.start()]
    rest = segment[match.end() :]
    token_match = _DRUG_TOKEN_AFTER_DOSAGE.match(rest)
    if not token_match:
        return _collect_drugs_from_segment(prefix)
    drug_token = token_match.group(1).strip()
    tail = rest[token_match.end() :]
    drugs: list[str] = []
    if prefix.strip() and not _is_duration_only_prefix(prefix):
        drugs.extend(_collect_drugs_from_segment(prefix))
    if drug_token:
        drugs.append(_label_part_drug(drug_token))
    drugs.extend(_collect_drugs_from_segment(tail.lstrip()))
    return drugs


def extract_drug_from_treatment(treatment: str | None) -> str | None:
    """Derive drug name from a parsed treatment string."""
    if treatment is None or not str(treatment).strip():
        return None
    text = str(treatment).strip()
    if " " not in text:
        return text
    text = _strip_leading_duration(text)
    name_part_comma = _strip_leading_duration(text.split(",", 1)[0].strip())
    parts: list[str] = []
    match = _DOSAGE_PREFIX.match(text)
    if match:
        name = match.group(1).strip()
        if "," in name:
            name = name.split(",", 1)[0].strip()
        parts = _collect_drugs_from_segment(name)
    elif _DOSAGE_INLINE.search(name_part_comma):
        parts = _collect_drugs_from_segment(name_part_comma)
    for label in _inhibitor_drug_labels(name_part_comma):
        parts = _merge_drug_parts(parts, label)
    return ", ".join(parts) if parts else None


def create_schema(con) -> None:
    con.execute(_DDL)


def upsert_diagnosis(con, name: str) -> int:
    """Return the diagnosis_id for `name`, creating the row if needed."""
    normalized = name.strip().lower()
    if not normalized:
        raise ValueError("diagnosis name must not be empty")
    con.execute(
        "INSERT INTO diagnosis (diagnosis_id, diagnosis_name) "
        "VALUES (nextval('diagnosis_id_seq'), ?) ON CONFLICT DO NOTHING",
        [normalized],
    )
    return con.execute(
        "SELECT diagnosis_id FROM diagnosis WHERE diagnosis_name = ?", [normalized]
    ).fetchone()[0]


def open_matrix(path: str) -> io.TextIOBase:
    opener = gzip.open if path.endswith(".gz") else open
    return opener(path, "rt", encoding="utf-8", errors="replace")


def read_matrix_header(path: str) -> tuple[dict[str, list[str]], list[list[str]]]:
    """Read `!` header fields and every Sample_characteristics_ch1 row."""
    fields: dict[str, list[str]] = {}
    characteristics_rows: list[list[str]] = []
    with open_matrix(path) as handle:
        for line in handle:
            if line.startswith(TABLE_BEGIN_MARKER):
                break
            if not line.startswith("!"):
                continue
            key, *values = next(csv.reader([line], delimiter="\t", quotechar='"'))
            key = key.lstrip("!")
            if key == "Sample_characteristics_ch1":
                characteristics_rows.append(values)
            elif key in KEY_COLUMNS or key in WANTED_KEYS:
                fields[key] = values
    return fields, characteristics_rows


def read_header(path: str) -> dict[str, list[str]]:
    """Read only scalar header keys (not characteristics rows)."""
    fields, _ = read_matrix_header(path)
    return fields


def missing_keys(fields: dict[str, list[str]], strict: bool = False) -> list[str]:
    required = (
        WANTED_KEYS if strict else tuple(k for k in WANTED_KEYS if k not in OPTIONAL_KEYS)
    )
    return [key for key in required if not fields.get(key)]


def _scalar(fields: dict[str, list[str]], key: str) -> str | None:
    values = fields.get(key)
    if not values or not str(values[0]).strip():
        return None
    return values[0]


def _as_int(value) -> int | None:
    if value is None:
        return None
    match = re.search(r"-?\d+", str(value))
    return int(match.group()) if match else None


def _cast(column: str, value):
    return _as_int(value) if column in INTEGER_COLUMNS else value


def parse_matrix(path: str, diagnosis_id: int, strict: bool = False):
    """Return (dataset_row, sample_rows) for one series matrix file.

    Raises MatrixError when required keys are missing so the caller can report
    the file instead of writing a partial study.
    """
    fields, characteristics_rows = read_matrix_header(path)
    absent = missing_keys(fields, strict=strict)
    if absent:
        raise MatrixError("missing required keys: " + ", ".join(absent))

    filename_platform = re.search(r"(?:^|-)GPL([0-9]+)(?:_series_matrix|$)", os.path.basename(path))
    series_platform = (
        f"GPL{filename_platform.group(1)}"
        if filename_platform
        else _scalar(fields, "Series_platform_id")
    )
    dataset_row = (
        diagnosis_id,
        os.path.basename(path),
        *(_scalar(fields, key) for key in SERIES_KEYS),
    )
    platform_index = SERIES_KEYS.index("Series_platform_id")
    dataset_row = dataset_row[:2 + platform_index] + (series_platform,) + dataset_row[3 + platform_index:]

    n_samples = len(fields["Sample_geo_accession"])

    def per_sample(key: str) -> list:
        values = fields.get(key) or []
        if len(values) < n_samples:
            values = values + [None] * (n_samples - len(values))
        return values[:n_samples]

    header_columns = {key: per_sample(key) for key in HEADER_SAMPLE_KEYS}
    treatments = extract_treatment_from_characteristics(characteristics_rows, n_samples)
    controls = extract_control_from_characteristics(characteristics_rows, n_samples)

    sample_rows = []
    for i in range(n_samples):
        header_part = tuple(
            _cast(KEY_COLUMNS[key], header_columns[key][i]) for key in HEADER_SAMPLE_KEYS
        )
        treatment = treatments[i]
        drug = extract_drug_from_treatment(treatment)
        sample_rows.append(header_part + (treatment, controls[i], drug))
    return dataset_row, sample_rows


def series_accession(path: str) -> str | None:
    return _scalar(read_header(path), "Series_geo_accession")
