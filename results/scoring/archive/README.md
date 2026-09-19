# Archived figures

Superseded outputs, kept for provenance. None of these should be shared or presented:
each has been replaced by something built on more data or better controls.

**Every file here carries a burned-in `SUPERSEDED` watermark and a footer naming its
replacement.** A README does not travel with a PDF that gets dragged into a slide deck, and
these figures look as finished as the current ones, so the warning has to live inside the
file. Re-apply with `tools/stamp_archive.sh`; unstamped copies stay in `originals/`, which
is gitignored.

| File | Why archived | Replaced by |
| --- | --- | --- |
| `figure_drug_ranking.pdf` / `.png` | Forest plot built on the 7-study pilot (20 contrasts). Sound design, but it shows a subset of the evidence and reports raw P values without the permutation result. | `../figure_singscore.pdf` (254 contrasts, rank-based) |
| `drug_ranking_by_cellline.png` | Bar chart of the same pilot contrasts with no confidence intervals. Bar length reads as confidence, which it is not. | `../figure_singscore.pdf` |
| `drug_ranking.png` | The first ranking, collapsed across cell lines. Averaging hid the two LSD1 arms cancelling to zero, which was the point of splitting. | `../figure_singscore.pdf` |
| `target_aml_km_by_age.png` | Day 1 baseline, survival by age band. The age-adjustment analysis states the same point more rigorously (age C-index 0.550, P = 0.05, and pLSC6 independent of age at P = 1.2e-5). | `R/incremental_value.R` output |

| `figure_combined_ranking.pdf` | The 18 weighted-sum contrasts that beat a random signature, with a scramble shRNA negative control ranked first. Superseded once rank-based enrichment showed the weighted sum was the problem, not the gene sets. **Worth retrieving if anyone asks how we knew the weighted-sum tail was noise** &mdash; it shows that in one image. | `../figure_singscore.pdf` |
| `figure_corpus_ranking.pdf` | Weighted-sum forest plot, array corpus only. | `../figure_singscore.pdf` |

## Pilot-era result tables — `pilot/`

The first working pass covered 4 to 6 studies and 20 contrasts. Those tables were named
`sample_scores.csv`, `drug_ranking.csv`, `contrast_statistics.csv`, `null_signature_test.csv`
and `null_drug_test.rds`, and sat in `results/` next to the corpus-scale files, where nothing
in the filename said which was which. `drug_ranking.csv` in particular reads like the main
deliverable and is not. They are here so the distinction is structural rather than a matter of
remembering. The corpus-scale equivalents are `corpus_drug_ranking.csv` (95 contrasts),
`ncbi_drug_ranking.csv` (162) and `lsc18_singscore_ranking.csv` (254, the primary result).

## Currently shareable

| File | Role |
| --- | --- |
| `../figure_null_tests.pdf` | **Primary.** Three permutation nulls: the score is specific in patients, and a weighted sum is not specific in either cell-line corpus. This motivates the change of summarisation shown in `figure_singscore.pdf`. |
| `../plsc6_target_km.png` | **Primary.** Survival curves separating by pLSC6 group. |
| `../figure_singscore.pdf` | **Primary.** Rank-based enrichment recovers signal the weighted sum could not, and the surviving agents group by mechanism. |

The three primary figures carry the whole argument. The forest plot is for anyone who asks
for per-contrast numbers.
