# Figure captions

Journal-ready captions for the three shareable figures, in narrative order: the score is
prognostic in patients (Fig 1), a weighted sum of it carries no specific signal in cell lines
(Fig 2), and rank-based enrichment of the gene set recovers that signal and identifies agents
that group by mechanism (Fig 3). Archived figures in `results/archive/` are superseded and
carry no caption.

Captions state what is shown and how it was computed. Interpretation belongs in the text.

---

## Figure 1 — `plsc6_target_km.pdf`

> **Figure 1. Overall survival in TARGET-AML stratified by the pLSC6 leukemic stem cell
> score.** Kaplan-Meier estimates of overall survival for 351 patients from the TARGET-AML
> cohort (NCI Genomic Data Commons, open-access tier) with both diagnostic RNA-seq and
> follow-up data; 134 deaths were observed. Each patient was scored with the six-gene
> pediatric leukemic stem cell signature (pLSC6; Elsayed et al., 2020) applied to
> log<sub>2</sub>(FPKM + 1) expression, and the cohort was dichotomised at the 60th
> percentile to match the split reported in the source publication (low, *n* = 211; high,
> *n* = 140). Five-year overall survival was 65.9% (95% CI 59.4–73.2) in the low-score group
> and 46.8% (38.6–56.9) in the high-score group. Cox proportional-hazards regression gave a
> hazard ratio of 1.83 (95% CI 1.30–2.56, *P* = 5.2 × 10<sup>−4</sup>) for the high versus
> low group and 1.90 (1.43–2.53, *P* = 9.3 × 10<sup>−6</sup>) per unit of the continuous
> score; Harrell's concordance index was 0.602. Tick marks denote censored observations.

*Optional sentence if the reviewer needs the comparison to the source publication:* The
published validation in an overlapping TARGET cohort (*n* = 205) reported a hazard ratio of
2.81 (1.85–4.28); the attenuation here is consistent with the use of a fixed percentile
rather than a recursive-partitioning cut point and a broader case selection.

---

## Figure 2 — `figure_null_tests.pdf`

> **Figure 2. A weighted sum of the signature is prognostic in patients but carries no
> specific signal as a pharmacodynamic readout in cell lines.** Each panel compares the
> observed statistic (vertical line, arrowhead) against a null distribution of 1,000 random
> six-gene signatures carrying the published pLSC6 coefficients (histogram). Shading marks
> the tail at least as extreme as the observed value; empirical *P* values use the
> (1 + *k*)/(*N* + 1) estimator.
> **(A)** Survival in 351 TARGET-AML patients, Cox *z* statistic: *z* = 4.43 against a null
> median of 1.24, *P* = 0.005. Note that 31% of random six-gene signatures were themselves
> significantly associated with survival at *P* < 0.05.
> **(B)** Drug response across 95 contrasts from 38 array-based series: 0.576 observed
> against a null median of 0.692, *P* = 0.99.
> **(C)** The same statistic across 162 contrasts from 38 further series quantified as
> NCBI-recomputed RNA-seq counts: 0.677 against 0.654, *P* = 0.39.
> Panels B and C use a weighted sum of expression, which carries no internal reference and is
> therefore moved by any broad transcriptional perturbation. They motivate the change of
> summarisation shown in Figure 3, and should be read as the diagnosis rather than the
> conclusion. Green marks an observed value outside its null; red marks one inside.

## Figure 3 — `figure_singscore.pdf` *(primary figure)*

> **Figure 3. Rank-based enrichment of the leukemic stem cell gene set recovers signal that a
> weighted sum does not, and the agents identified group by mechanism.**
> **(A)** Bar height gives the number of treated-versus-control contrasts surviving
> Benjamini-Hochberg correction on the permutation *P*, across 254 contrasts from 75 GEO
> series; the multiplier above each bar is the ratio of contrasts exceeding their null at
> raw *P* < 0.05 to the number expected by chance at that threshold. The data, the contrasts
> and the permutation design are identical in every bar; only the summarisation and the gene
> set change. A weighted sum of pLSC6 yields no survivors and 1.4-fold enrichment. Rank-based
> enrichment (singscore; Foroutan et al., 2018), which ranks each sample against its own
> transcriptome so that a global shift cancels, yields 6 survivors at 3.0-fold for the same
> six genes, 23 at 4.2-fold for the 17-gene LSC17 set, and 31 at 4.2-fold for their 18-gene
> union. Because a rank-based score requires set membership and not coefficients, the
> LSC17 gene list is usable here although its published coefficients could not be verified.
> Across 2,000 permutations no random set of matched size reached the observed mean absolute
> effect for any of the three sets (*P* < 0.0005), and the result holds separately in each
> corpus: 23 of 92 array-based contrasts exceeded their null against 4.6 expected, and 30 of
> 162 counts-based contrasts against 8.1 expected.
> **(B)** All 20 contrasts that both lowered the LSC18 enrichment score and survived
> correction, grouped by agent, with Benjamini-Hochberg *q* shown per contrast and colour
> denoting the source of the expression data. Effects are in singscore units, bounded at
> ±0.5. The menin-MLL inhibitor BAY-1251152 lowered the score in seven independent cell lines,
> four of which carry *KMT2A* rearrangements, the dependency the compound is designed against.
> The largest single effect was the IDH1 inhibitor BAY1436032, which also survived in
> combination with azacitidine. Eleven further contrasts survived correction in the opposite
> direction and are not shown.

## Notes on use

**Abbreviations to define at first use in the main text:** pLSC6, six-gene pediatric
leukemic stem cell score; GEO, Gene Expression Omnibus; FPKM, fragments per kilobase of
transcript per million mapped reads; CI, confidence interval; FDR, false discovery rate.

**Gene nomenclature.** *KMT2A* is the current HGNC symbol for the gene the source studies
annotate as MLL. Figure 2C and Figure 3 reproduce each study's own label; use *KMT2A*
(MLL) in text. Similarly *ADGRG1* is the current symbol for GPR56, which appears in the
published pLSC6 equation and in `data/signatures/plsc6_elsayed2020.csv`.

**Statistics.** Permutation *P* values in Figure 2 use 1,000 draws and those in Figure 3 use
2,000, with the (1 + *k*)/(*N* + 1) estimator, so the smallest reportable values are 0.001 and
0.0005 respectively. Figure 3's global tests reached that floor: no permutation among 2,000
matched the observed value for any gene set. No adjustment for multiple comparisons
was applied to the permutation *P* values shown in Figure 3; the single contrast surviving
that adjustment is named in the caption.

**What not to claim from these figures.** On cell lines the score measures a
stemness-associated expression programme, not prognosis; prognostic language applies only to
Figure 1. Effects in Figure 3B are small in absolute terms, most between −0.04 and −0.09 on a
scale bounded at ±0.5, and eleven surviving contrasts move in the opposite direction. The
LSC17 coefficients remain unverified; only its gene membership is used. And the pinometostat
result that survived under the weighted sum **did not replicate** on independent data and is
not among the agents reported here.
