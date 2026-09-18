"""Load GEO sample metadata into a DuckDB file for a diagnosis's downloaded matrices.

Intended to run against the `matrices/` directory produced by
R_Scripts/05_Download_Metadata_Inventory.R for a given diagnosis, e.g.:
    python parsing/buildDb.py downloads/geo_aml/matrices --db-path data/geo.duckdb

When downloaded_matrices.csv (or all_results.csv) exists next to data_dir, it is
loaded into a `studies` table (study-level inventory from run_geo_pipeline()).
The `samples` table is always built by parsing *_series_matrix.txt.gz files under
data_dir for per-sample fields (series_accession, sample_geo_accession, etc.).

Multi-platform studies ship one matrix per platform (`GSE*-GPL*_series_matrix.txt.gz`).
`series_platform_id` comes from the GPL in that filename (same rule as R/bronze.R
`platform_from_file`), not from duplicate `!Series_platform_id` header lines.

After changing this module, rebuild the Shiny database from cached matrices, e.g.:
    python parsing/buildDb.py downloads/geo_aml/matrices --db-path data/geo.duckdb
"""

import argparse
import csv
import glob
import gzip
import os
import re

import duckdb

REPO_ROOT = os.path.join(os.path.dirname(__file__), "..")
DEFAULT_DB_PATH = os.path.join(REPO_ROOT, "data", "geo.duckdb")
TABLE_BEGIN_MARKER = "!series_matrix_table_begin"
REPORT_FILENAMES = ("downloaded_matrices.csv", "all_results.csv")

WANTED_KEYS = (
    "Series_geo_accession",
    "Series_platform_id",
    "Sample_geo_accession",
    "Sample_organism_ch1",
    "Sample_data_row_count",
    "Sample_molecule_ch1",
    "Series_pubmed_id",
)


def parse_tsv_line(line):
    return next(csv.reader([line], delimiter="\t", quotechar='"'))


def platform_from_matrix_path(path):
    """GPL from matrix basename, e.g. GSE100708-GPL16791_series_matrix.txt.gz -> GPL16791."""
    match = re.search(r"GPL\d+", os.path.basename(path))
    return match.group(0) if match else None


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


_EXCLUDED_TREATMENT_VALUES = frozenset({"untreated", "vehicle", "dmso", "control"})


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
    """Parsed value from a treatment/drug characteristics cell (before control vs drug split)."""
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
    """DMSO, vehicle, untreated, and control-like values from the same characteristics row."""
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


def parse_series_matrix(path):
    fields = {}
    characteristics_rows = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(TABLE_BEGIN_MARKER):
                break  # stop before the (large) expression matrix
            if not line.strip().startswith("!"):
                continue
            key, *values = parse_tsv_line(line)
            key = key.lstrip("!")
            if key == "Sample_characteristics_ch1":
                characteristics_rows.append(values)
            elif key in WANTED_KEYS:
                fields[key] = values

    series_accession = fields["Series_geo_accession"][0]
    platform_id = platform_from_matrix_path(path) or fields["Series_platform_id"][0]
    pubmed_raw = fields.get("Series_pubmed_id")
    series_pubmed_id = (
        pubmed_raw[0] if pubmed_raw and str(pubmed_raw[0]).strip() else None
    )
    sample_ids = fields["Sample_geo_accession"]
    n_samples = len(sample_ids)

    def per_sample(key):
        values = fields.get(key)
        if not values:
            return [None] * n_samples
        if len(values) < n_samples:
            values = values + [None] * (n_samples - len(values))
        return values[:n_samples]

    organisms = per_sample("Sample_organism_ch1")
    row_counts = per_sample("Sample_data_row_count")
    molecules = per_sample("Sample_molecule_ch1")
    treatments = extract_treatment_from_characteristics(characteristics_rows, n_samples)
    controls = extract_control_from_characteristics(characteristics_rows, n_samples)

    rows = []
    for sample_id, organism, row_count, treatment, control, molecule in zip(
        sample_ids, organisms, row_counts, treatments, controls, molecules
    ):
        parsed_count = None
        if row_count is not None and str(row_count).strip():
            try:
                parsed_count = int(row_count)
            except ValueError:
                parsed_count = None
        rows.append(
            (
                series_accession,
                platform_id,
                series_pubmed_id,
                sample_id,
                organism,
                parsed_count,
                treatment,
                control,
                molecule,
            )
        )
    return rows


def find_report_path(data_dir):
    """Find run_geo_pipeline()'s inventory report next to a matrices/ directory."""
    parent = os.path.dirname(os.path.normpath(data_dir))
    for name in REPORT_FILENAMES:
        candidate = os.path.join(parent, name)
        if os.path.exists(candidate):
            return candidate
    return None


def load_samples_from_matrices(con, data_dir):
    con.execute(
        "CREATE TABLE samples ("
        "series_accession VARCHAR, series_platform_id VARCHAR, "
        "series_pubmed_id VARCHAR, "
        "sample_geo_accession VARCHAR, sample_organism_ch1 VARCHAR, "
        "sample_data_row_count INTEGER, treatment VARCHAR, control VARCHAR, "
        "sample_molecule_ch1 VARCHAR)"
    )

    paths = sorted(glob.glob(os.path.join(data_dir, "*_series_matrix.txt.gz")))
    for path in paths:
        rows = parse_series_matrix(path)
        con.executemany("INSERT INTO samples VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", rows)
        print(f"Loaded {rows[0][0]}: {len(rows)} samples")

    return len(paths)


def build_db(data_dir, db_path=None):
    """Build studies (optional) and samples tables for the matrices in data_dir.

    db_path=None keeps the database in memory; otherwise the file at db_path is
    deleted first so every build starts from a clean database.
    """
    if db_path:
        os.makedirs(os.path.dirname(db_path), exist_ok=True)
        if os.path.exists(db_path):
            os.remove(db_path)
        con = duckdb.connect(db_path)
    else:
        con = duckdb.connect()

    report_path = find_report_path(data_dir)
    if report_path:
        con.execute("CREATE TABLE studies AS SELECT * FROM read_csv_auto(?)", [report_path])
        count = con.execute("SELECT COUNT(*) FROM studies").fetchone()[0]
        print(f"Loaded studies from {report_path}: {count} rows")

    matrix_count = load_samples_from_matrices(con, data_dir)
    sample_count = con.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    print(f"Parsed {matrix_count} matrix file(s); {sample_count} sample row(s) in samples")

    return con


def main():
    parser = argparse.ArgumentParser(
        description="Build the GEO samples DuckDB from series matrix files."
    )
    parser.add_argument(
        "data_dir",
        help="Directory of *_series_matrix.txt.gz files, e.g. downloads/geo_aml/matrices.",
    )
    parser.add_argument(
        "--db-path", default=DEFAULT_DB_PATH,
        help="DuckDB file to write; deleted and rebuilt from scratch each run.",
    )
    args = parser.parse_args()

    con = build_db(args.data_dir, args.db_path)
    print(con.execute("SELECT COUNT(*) FROM samples").fetchone())
    con.close()


if __name__ == "__main__":
    main()

