#!/usr/bin/env Rscript
# Rank-based single-sample enrichment scoring (singscore, Foroutan et al. 2018).
#
# Why this exists. The published scores are weighted sums of expression, which carry no
# correction for the rest of the transcriptome. A drug that perturbs expression broadly
# moves any gene set, which is exactly what the permutation tests detected: random six-gene
# sets separated treated from control samples as well as pLSC6 did.
#
# A rank-based score asks a different question: are the set's genes unusually HIGH relative
# to everything else measured in that same sample? Ranking is computed within each sample,
# so a global shift cancels, and the score is bounded and comparable across samples without
# further standardisation.
#
# It also removes the coefficient problem. An enrichment score needs the gene SET, not the
# weights, so the LSC17 gene list is usable even though its published coefficients could not
# be verified.

# ---- gene sets -------------------------------------------------------------

# LSC17 (Ng et al. 2016) gene membership. The COEFFICIENTS remain unverified and are not
# used here; only the identity of the genes matters for a rank-based score.
LSC17_GENES <- c("DNMT3B","ZBTB46","NYNRIN","ARHGAP22","LAPTM4B","MMRN1","DPYSL3",
                 "KIAA0125","CDK6","CPXM1","SOCS2","SMIM24","EMP1","NGFRAP1",
                 "CD34","AKR1C3","GPR56")
PLSC6_GENES <- c("DNMT3B","GPR56","CD34","SOCS2","SPINK2","FAM30A")

# KIAA0125 and FAM30A are the same locus, so the union is 18 genes.
LSC18_GENES <- unique(c(LSC17_GENES, setdiff(PLSC6_GENES, c("FAM30A"))))


# Entrez ids for every LSC gene, resolved against NCBI Gene rather than assumed. Needed
# because NCBI-recomputed count matrices are keyed on GeneID, not symbol; without this the
# 17- and 18-gene sets silently lose every sequencing study to a coverage threshold.
LSC_ENTREZ <- c(
  AKR1C3 = "8644",
  ARHGAP22 = "58504",
  CD34 = "947",
  CDK6 = "1021",
  CPXM1 = "56265",
  DNMT3B = "1789",
  DPYSL3 = "1809",
  EMP1 = "2012",
  FAM30A = "9834",
  GPR56 = "9289",
  KIAA0125 = "9834",
  LAPTM4B = "55353",
  MMRN1 = "22915",
  NGFRAP1 = "27018",
  NYNRIN = "57523",
  SMIM24 = "284422",
  SOCS2 = "8835",
  SPINK2 = "6691",
  ZBTB46 = "140685"
)

#' Relabel a counts matrix keyed on Entrez ids so LSC gene symbols resolve.
relabel_entrez <- function(M) {
  inv <- setNames(names(LSC_ENTREZ), unname(LSC_ENTREZ))
  hit <- which(rownames(M) %in% names(inv))
  if (length(hit)) rownames(M)[hit] <- inv[rownames(M)[hit]]
  M
}

GENE_SETS <- list(pLSC6 = PLSC6_GENES, LSC17 = LSC17_GENES, LSC18 = LSC18_GENES)

# Aliases seen across GEO platforms and the NCBI counts annotation.
GENE_ALIASES <- list(
  GPR56    = c("ADGRG1"), FAM30A = c("KIAA0125"), KIAA0125 = c("FAM30A"),
  NGFRAP1  = c("BEX3"),   SMIM24 = c("C19orf77","HSPC323"), ZBTB46 = c("BTBD4"),
  NYNRIN   = c("KIAA1305","CGIN1"), CPXM1 = c("CPXM"), DPYSL3 = c("CRMP4")
)

resolve_set <- function(genes, available) {
  vapply(genes, function(g) {
    if (g %in% available) return(g)
    for (a in GENE_ALIASES[[g]]) if (!is.null(a) && a %in% available) return(a)
    NA_character_
  }, character(1))
}

# ---- the score -------------------------------------------------------------

#' Rank a matrix once, within sample. High expression takes a high rank.
#' Ranking dominates the cost of singscore and is independent of the gene set, so it is
#' computed once per study and reused across the observed score and every permutation.
sing_ranks <- function(expr) apply(expr, 2, rank, ties.method = "average")

#' singscore for an up-regulated gene set.
#'
#' Ranks each sample independently, takes the mean rank of the set, and rescales against
#' the minimum and maximum mean rank the set could attain in a transcriptome of that size.
#' Returns values on [-0.5, 0.5]; zero is the rank a random set of the same size expects.
#'
#' @param expr   genes x samples, any monotone transform of expression
#' @param set    gene symbols; those absent are dropped and reported
#' @param min_frac  minimum fraction of the set that must be present
#' @param R      optional precomputed rank matrix from sing_ranks(); ranking is by far the
#'               expensive step and does not depend on the gene set, so a permutation loop
#'               should rank once and pass it in rather than re-ranking every draw.
sing_score <- function(expr, set, min_frac = 0.8, R = NULL) {
  present <- resolve_set(set, rownames(expr))
  found   <- present[!is.na(present)]
  if (length(found) / length(set) < min_frac)
    stop(sprintf("Only %d/%d set genes present (%.0f%%).",
                 length(found), length(set), 100*length(found)/length(set)), call. = FALSE)

  n <- nrow(expr); m <- length(found)
  if (is.null(R)) R <- sing_ranks(expr)
  obs <- colMeans(R[found, , drop = FALSE])

  lo <- mean(seq_len(m))                 # set occupies the lowest ranks
  hi <- mean(seq(n - m + 1, n))          # set occupies the highest ranks
  structure((obs - lo) / (hi - lo) - 0.5,
            genes_found = m, genes_expected = length(set), symbols = found)
}

#' Dispersion of the set's ranks, which singscore reports alongside the score.
#' A low score with high dispersion means the set is not behaving coherently.
sing_dispersion <- function(expr, set, min_frac = 0.8) {
  present <- resolve_set(set, rownames(expr)); found <- present[!is.na(present)]
  if (!length(found)) return(rep(NA_real_, ncol(expr)))
  R <- apply(expr, 2, function(v) rank(v, ties.method = "average"))
  apply(R[found, , drop = FALSE], 2, stats::mad)
}
