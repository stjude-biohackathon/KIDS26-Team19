"""Migrate an older geo.db that incorrectly made GSE accession unique.

    python parsing/migrateDb.py --db-path data/geo.db

The migration preserves diagnosis, dataset, and sample rows while replacing the
dataset table constraint with source-file uniqueness, allowing multi-GPL studies.
Sample rows gain treatment/control derived from legacy sample_characteristics_ch1 when present.
"""
from __future__ import annotations

import argparse
import os
import tempfile

import duckdb

from schema import create_schema, derive_treatment_control_from_cell


def migrate(db_path: str) -> None:
    db_path = os.path.abspath(db_path)
    if not os.path.exists(db_path):
        raise FileNotFoundError(db_path)
    parent = os.path.dirname(db_path)
    fd, replacement = tempfile.mkstemp(prefix="geo-migration-", suffix=".db", dir=parent)
    os.close(fd)
    os.remove(replacement)
    try:
        source = duckdb.connect(db_path, read_only=True)
        target = duckdb.connect(replacement)
        escaped = db_path.replace("'", "''")
        target.execute(f"ATTACH '{escaped}' AS old (READ_ONLY)")
        create_schema(target)
        target.execute("INSERT INTO diagnosis (diagnosis_name) SELECT diagnosis_name FROM old.diagnosis")
        target.execute("""
            INSERT INTO dataset (diagnosis_id, source_file, series_geo_accession,
                                 series_platform_id, series_pubmed_id)
            SELECT nd.diagnosis_id, od.source_file, od.series_geo_accession,
                   od.series_platform_id, od.series_pubmed_id
            FROM old.dataset od
            JOIN old.diagnosis odg ON odg.diagnosis_id = od.diagnosis_id
            JOIN diagnosis nd ON nd.diagnosis_name = odg.diagnosis_name
        """)
        old_sample_cols = {row[0] for row in target.execute("DESCRIBE old.sample").fetchall()}
        has_legacy_chars = "sample_characteristics_ch1" in old_sample_cols
        char_select = (
            "os.sample_characteristics_ch1"
            if has_legacy_chars
            else "CAST(NULL AS VARCHAR)"
        )
        sample_rows = target.execute(f"""
            SELECT nd.dataset_id, nd.diagnosis_id, os.sample_geo_accession,
                   os.sample_organism_ch1, os.sample_data_row_count,
                   {char_select}, os.sample_molecule_ch1
            FROM old.sample os
            JOIN old.dataset od ON od.dataset_id = os.dataset_id
            JOIN dataset nd ON nd.source_file = od.source_file
        """).fetchall()
        for dataset_id, diag_id, gsm, org, row_count, chars, molecule in sample_rows:
            treatment, control = derive_treatment_control_from_cell(chars)
            target.execute(
                """
                INSERT INTO sample (
                    dataset_id, diagnosis_id, sample_geo_accession,
                    sample_organism_ch1, sample_data_row_count,
                    treatment, control, sample_molecule_ch1
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [dataset_id, diag_id, gsm, org, row_count, treatment, control, molecule],
            )
        target.close()
        source.close()
        os.replace(replacement, db_path)
    except Exception:
        if os.path.exists(replacement):
            os.remove(replacement)
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--db-path", default="data/geo.db")
    args = parser.parse_args()
    migrate(args.db_path)
    print(f"Migrated {args.db_path}; multiple platforms per GSE are now allowed.")
