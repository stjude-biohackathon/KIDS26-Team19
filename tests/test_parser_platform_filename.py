import gzip
from pathlib import Path

from parsing.schema import parse_matrix


def test_filename_platform_overrides_header(tmp_path: Path):
    path = tmp_path / "GSE14-GPL11_series_matrix.txt.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(
            '!Series_geo_accession\t"GSE14"\n'
            '!Series_platform_id\t"GPL1485"\n'
            '!Sample_geo_accession\t"GSM1"\n'
            '!Sample_organism_ch1\t"Homo sapiens"\n'
            '!Sample_data_row_count\t"100"\n'
            '!Sample_molecule_ch1\t"RNA"\n'
            '!series_matrix_table_begin\n'
        )
    dataset, samples = parse_matrix(str(path), diagnosis_id=1)
    assert dataset[3] == "GPL11"
    assert len(samples) == 1
