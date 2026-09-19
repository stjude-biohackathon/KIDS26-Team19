# Leukemic stem cell signatures

Two scores live here. Both are stored as editable CSVs so the numbers can be checked
against the primary source rather than buried in code.

| File | Score | Source |
| --- | --- | --- |
| `lsc17_ng2016.csv` | LSC17, 17 genes | Ng SWK et al. *A 17-gene stemness score for rapid determination of risk in acute leukaemia.* Nature 2016;540:433-437. doi:10.1038/nature20598 |
| `plsc6_elsayed2020.csv` | pLSC6, 6 genes, pediatric | Elsayed AH et al. *A six-gene leukemic stem cell score identifies high risk pediatric acute myeloid leukemia.* Leukemia 2020;34:735-745. doi:10.1038/s41375-019-0604-8 |

## Verification status

| Signature | Status | Checked against |
| --- | --- | --- |
| **pLSC6** | **VERIFIED 2026-09-16** | Elsayed et al., Leukemia 2020, full text via PubMed Central (PMC7135934) |
| **LSC17** | **NOT VERIFIED** | Ng et al., Nature 2016 is paywalled and no open source reproduces the coefficient table |

### pLSC6 was wrong and is now fixed

The first draft of `plsc6_elsayed2020.csv` had all six coefficients attached to the
**wrong genes** -- the correct six values, permuted. The paper states the equation
verbatim:

> pLSC6 = (DNMT3B x 0.189) + (GPR56 x 0.054) + (CD34 x 0.0171) + (SOCS2 x 0.141)
> + (SPINK2 x 0.109) + (FAM30A x 0.0516)

| Gene | Was (wrong) | Now (published) |
| --- | --- | --- |
| DNMT3B | 0.0171 | **0.189** |
| GPR56 | 0.109 | **0.054** |
| CD34 | 0.141 | **0.0171** |
| SOCS2 | 0.0516 | **0.141** |
| SPINK2 | 0.054 | **0.109** |
| FAM30A | 0.189 | **0.0516** |

Every score computed before this correction was wrong, and nothing in the output would
have indicated it. This is the reason the `verified` flag and the refusal in
`score_lsc()` exist, and the reason the LSC17 **coefficients** must not be used until
someone repeats this exercise against the Nature paper.

**This restricts the weighted sum only.** A rank-based (singscore) enrichment uses gene
membership and never touches the coefficients, so LSC17 and LSC18 are safe to use there and
every reported LSC17/LSC18 result comes from that path. `score_lsc(sig = "LSC17")` still
errors by design; `sing_score()` does not.

### LSC17 still needs a human

`lsc17_ng2016.csv` is unverified and its coefficients came from the same unreliable
source that got pLSC6 wrong, so **assume they are wrong until checked**. St. Jude will
have institutional access to Nature 2016;540:433-437, and Dr. Stan Pounds is a
co-author on the pLSC6 paper and is consulting on this project. Either route is faster
and more reliable than anything reachable from here.

## Scale requirement

Both scores were derived on **log-scale** expression. Elsayed et al. used log-transformed
MAS 5.0 U133A signals for discovery and `log2(RPKM + 1)` for the TARGET validation.
Feeding an un-logged matrix produces values on a completely different scale.
`score_lsc()` flags this and reports `InputScale` per row.

## What these scores mean, and what they do not

Both scores were derived from **patient** samples and validated as **prognostic**:
they stratify overall survival in AML patients. That is the only setting in which
the word "prognostic" applies to them.

**On a cell line, neither score is prognostic.** A cell line has no survival outcome.
Scoring drug-treated cell lines measures a *shift in a stemness-associated expression
signature*, which is a pharmacodynamic readout. That is a legitimate and interesting
thing to measure -- "does this drug move cells away from an LSC-like state?" -- but it
is a different claim from prognosis, and it needs different wording in the demo and
the README.

If the project wants to claim prognostic relevance, the chain has to be closed with
patient data: show that the score stratifies survival in a patient cohort, then show
that the drug moves the score in cell lines. TARGET-AML supplies the first half with
open-access expression and survival on the same subjects.

## Scale and comparability

The score is a weighted sum of expression values, so it inherits the units of the
input matrix. Scores are **not comparable across platforms, studies, or
normalizations**. Valid comparisons are within a single study, most usefully treated
versus vehicle on the same cell line and platform. `score_lsc()` returns a
within-study z-score alongside the raw score for this reason.
