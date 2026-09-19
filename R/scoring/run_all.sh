#!/usr/bin/env bash
# Rebuild every result from scratch. Expects R with `survival`, `data.table` and `testthat`.
# Network access is required: expression data is pulled from NCBI GEO and the NCI GDC.
# Scoring stages are resumable and skip studies already written to results/.
set -euo pipefail
cd "$(dirname "$0")"

echo "== unit tests =="
Rscript -e 'testthat::test_dir("tests/testthat")'

echo "== corpus composition (uses cached headers under bronze/headers) =="
Rscript R/scan_corpus.R

echo "== score the array route (series matrices) =="
Rscript R/score_corpus.R

echo "== score the RNA-seq route (NCBI-recomputed counts) =="
Rscript R/score_ncbi_corpus.R

echo "== contrasts: harmonise treatments, treated vs control =="
# The third argument is an output prefix, so corpus runs do not overwrite each other.
Rscript R/rank_drugs.R bronze/matrices pLSC6 corpus_   # -> results/corpus_drug_ranking.csv
Rscript R/rank_drugs.R Toy-Datasets    pLSC6           # -> results/drug_ranking.csv (worked example)

echo "== permutation tests, weighted sum (the diagnosis) =="
Rscript R/null_signature_survival.R 1000   # survival, 351 TARGET patients
Rscript R/null_corpus_test.R 1000          # 95 array contrasts
Rscript R/null_ncbi_test.R 1000            # 162 NCBI-counts contrasts

echo "== permutation test, rank-based enrichment (the recovery) =="
Rscript R/singscore_analysis.R 2000        # pLSC6, LSC17, LSC18 over 254 contrasts

echo "== survival models =="
Rscript R/validate_plsc6_target.R          # Figure 1, Kaplan-Meier
Rscript R/incremental_value.R              # nested Cox: does the score add beyond age

echo "== figures =="
Rscript R/figure_null_tests.R              # Figure 2, three permutation nulls
Rscript R/figure_singscore.R               # Figure 3, rank-based recovery

echo "done. Shareable figures and tables are in results/; superseded ones in results/archive/."
