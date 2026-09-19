# Scoring layer

Takes a GEO expression matrix, returns a per-sample stemness score, and tests whether that
score carries information the same number of random genes would not.

Sits downstream of the bronze layer: it reads the series matrices already in `bronze/`.
Run `./run_all.sh` for the whole thing in order, or the stages below individually.

## Read these first

| File | What it is |
| --- | --- |
| `lsc_scores.R` | **Library.** Signature loading, platform annotation, probe-to-gene collapsing, the weighted sum, and the refusal guard. Sourced by almost everything else. |
| `singscore.R` | **Library.** Rank-based enrichment. This is what the reported LSC17/LSC18 results use. |
| `harmonise.R` | **Library.** Pulls an agent and a control flag out of free-text sample characteristics. |
| `ncbi_counts.R` | **Library.** Reads the counts NCBI recomputes for GEO2R-flagged sequencing series. |

Those four define functions and are sourced. The rest are scripts you run.

## Pipeline order

| Stage | Script | Produces |
| --- | --- | --- |
| 1 | `scan_corpus.R` | `results/corpus_scan.csv` — which series support a treated/control design, from headers only |
| 2 | `score_corpus.R` | `results/corpus_sample_scores.csv` — array route, 38 studies, 613 samples |
| 3 | `score_ncbi_corpus.R` | `results/ncbi_sample_scores.csv` — RNA-seq route, 44 studies, 970 samples |
| 4 | `rank_drugs.R` | contrast tables (see the argument note below) |
| 5 | `null_signature_survival.R` | `results/null_survival_test.*` — permutation, survival |
| 5 | `null_corpus_test.R` | `results/corpus_null_test.csv` — permutation, 95 array contrasts |
| 5 | `null_ncbi_test.R` | `results/ncbi_null_test.csv` — permutation, 162 counts contrasts |
| 6 | `singscore_analysis.R` | `results/singscore_tests.rds` + `lsc18_singscore_ranking.csv` — the primary result |
| 7 | `validate_plsc6_target.R` | Figure 1, Kaplan-Meier on TARGET-AML |
| 7 | `incremental_value.R` | nested Cox models: does the score add beyond age |
| 8 | `figure_null_tests.R` | Figure 2 |
| 8 | `figure_singscore.R` | Figure 3 |

Stages 2, 3 and 6 are resumable: each study writes its own file and is skipped on a rerun.

## Things that will catch you out

**`rank_drugs.R` defaults to `Toy-Datasets`**, a 4-to-6 study pilot. Run bare, it writes
pilot-scale files. Pass a directory and an output prefix for corpus scale:

```bash
Rscript R/rank_drugs.R bronze/matrices pLSC6 corpus_   # -> results/corpus_drug_ranking.csv
```

**The LSC17 coefficients are unverified**, so `score_lsc(sig = "LSC17")` errors by design.
That restricts the weighted sum only. `sing_score()` needs gene membership, not weights, so
LSC17 and LSC18 are usable there, and every reported LSC17/LSC18 number comes from that path.
See `data/signatures/README.md`.

**`score_lsc()` refuses to run on a signature whose `verified` flag is FALSE.** That guard
exists because a draft coefficient table had all six pLSC6 values attached to the wrong genes,
and nothing in the output would have revealed it. Do not set the flag without checking the
primary paper.

**Scores are not comparable across studies in raw units.** The score inherits its input's
units, so every contrast is computed in within-study standard deviations.

**Counts matrices are keyed on Entrez GeneID, not gene symbols.** `singscore.R` carries the
mapping for all 19 genes. A partial relabel makes a set fall below its coverage threshold and
the study is dropped with no error.

## Dependencies

R 4.5.1 with `survival`, `data.table`, `testthat`, plus `curl` on the PATH. No Bioconductor,
deliberately. All inputs are open access: no credentials, no controlled-access data.

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'   # 26 assertions, no network
```
