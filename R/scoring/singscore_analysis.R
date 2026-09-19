#!/usr/bin/env Rscript
# Does rank-based enrichment recover signal that the weighted sum could not?
#
# Identical contrasts, identical permutation design; only the summarisation changes.
# Three gene sets are compared so that set size and summarisation are separable:
#   pLSC6  6 genes, the published prognostic set
#   LSC17 17 genes, membership only (its coefficients were never verifiable)
#   LSC18 18 genes, the union, the largest coherent LSC set available here
#
#   Rscript R/singscore_analysis.R [nperm]

set.seed(20260918)
source(file.path("R", "lsc_scores.R"))
source(file.path("R", "ncbi_counts.R"))
source(file.path("R", "singscore.R"))
NPERM <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) as.integer(a[1]) else 500L }

sig <- load_signature("pLSC6")
SETS <- GENE_SETS

# ---- assemble per-study matrices and contrast groups ------------------------

acc <- setNames(lapply(names(SETS), function(x)
        list(obs=c(), keys=c(), src=c(), null=list())), names(SETS))

process <- function(gse, M, meta, src) {
  d <- meta[match(colnames(M), meta$GSM), ]
  if (all(is.na(d$Agent))) return(invisible(NULL))
  d$Line <- ifelse(!is.na(d$CellLine), d$CellLine, gse)
  grp <- list()
  for (ln in unique(na.omit(d$Line))) {
    i  <- !is.na(d$Line) & d$Line == ln
    ct <- which(i & !is.na(d$Agent) & d$Agent == "control"); if (length(ct) < 2) next
    for (a in setdiff(unique(na.omit(d$Agent[i])), "control")) {
      tr <- which(i & !is.na(d$Agent) & d$Agent == a); if (length(tr) < 2) next
      grp[[length(grp)+1]] <- list(key = paste(gse, ln, a, sep="|"), tr = tr, ct = ct)
    }
  }
  if (!length(grp)) return(invisible(NULL))

  R <- sing_ranks(M)                       # ranked once, reused everywhere below
  pool <- which(apply(M, 1, stats::sd) > 0)
  ctr  <- function(sc) vapply(grp, function(g) mean(sc[g$tr]) - mean(sc[g$ct]), numeric(1))

  for (nm in names(SETS)) {
    o <- try(sing_score(M, SETS[[nm]], min_frac = 0.7, R = R), silent = TRUE)
    if (inherits(o, "try-error")) next
    m <- attr(o, "genes_found"); if (length(pool) < m * 5) next
    nl <- vapply(seq_len(NPERM), function(p)
            ctr(sing_score(M, rownames(M)[sample(pool, m)], min_frac = 0, R = R)),
          numeric(length(grp)))
    if (length(grp) == 1L) nl <- matrix(nl, nrow = 1)
    acc[[nm]]$obs  <<- c(acc[[nm]]$obs,  ctr(o))
    acc[[nm]]$keys <<- c(acc[[nm]]$keys, vapply(grp, function(g) g$key, character(1)))
    acc[[nm]]$src  <<- c(acc[[nm]]$src,  rep(src, length(grp)))
    acc[[nm]]$null[[length(acc[[nm]]$null)+1]] <<- t(nl)
  }
  message(sprintf("  %-11s %-6s %d genes, %d contrasts", gse, src, nrow(M), length(grp)))
  invisible(NULL)
}

nc <- read.csv(file.path("results","ncbi_sample_scores.csv"), stringsAsFactors = FALSE)
for (gse in unique(nc$GSE)) {
  M <- try(ncbi_counts_matrix(gse, sig, keep_all = TRUE)$expr, silent = TRUE)
  if (inherits(M, "try-error")) next
  M <- relabel_entrez(M)   # all 19 LSC symbols, not just the six with a signature row
  process(gse, M, nc[nc$GSE == gse, ], "counts")
  rm(M); gc(verbose = FALSE)
}

ar <- read.csv(file.path("results","corpus_sample_scores.csv"), stringsAsFactors = FALSE)
wanted <- unique(unlist(lapply(seq_len(nrow(sig)), function(i) symbols_for(sig[i,]))))
for (gse in unique(ar$GSE)) {
  mf <- file.path("bronze","matrices", paste0(gse,"_series_matrix.txt.gz"))
  if (!file.exists(mf)) next
  sm <- try(read_series_matrix(mf), silent = TRUE); if (inherits(sm,"try-error")) next
  if (isTRUE(max(sm$expr, na.rm=TRUE) > 50)) sm$expr <- log2(pmax(sm$expr,0)+1)
  pm <- try(fetch_platform_annotation(sm$gpl, prefer_symbols = wanted), silent = TRUE)
  if (inherits(pm,"try-error")) { rm(sm); gc(verbose=FALSE); next }
  pmk <- pm[!is.na(pm$probe) & !is.na(pm$symbol), ]
  pmk <- pmk[pmk$probe %in% rownames(sm$expr) & nzchar(pmk$symbol) &
             !grepl("///", pmk$symbol, fixed = TRUE), ]
  if (nrow(pmk) >= 500) {
    sub <- sm$expr[pmk$probe, , drop=FALSE]; mu <- rowMeans(sub, na.rm=TRUE)
    best <- unlist(tapply(seq_len(nrow(pmk)), pmk$symbol, function(i) i[which.max(mu[i])]))
    M <- sub[best,,drop=FALSE]; rownames(M) <- pmk$symbol[best]
    M <- M[stats::complete.cases(M), , drop=FALSE]
    process(gse, M, ar[ar$GSE == gse, ], "array")
    rm(M, sub)
  }
  rm(sm, pm, pmk); gc(verbose = FALSE)
}

run_set <- function(nm) {
  a <- acc[[nm]]; null <- do.call(cbind, a$null)
  emp <- vapply(seq_along(a$obs), function(j)
    (1 + sum(abs(null[,j]) >= abs(a$obs[j]))) / (NPERM + 1), numeric(1))
  list(set = nm, n = length(SETS[[nm]]), obs = a$obs, keys = a$keys, src = a$src,
       null_p = emp, global_obs = mean(abs(a$obs)), global_null = rowMeans(abs(null)))
}

res <- lapply(names(SETS), run_set); names(res) <- names(SETS)
saveRDS(res, file.path("results","singscore_tests.rds"))

# Per-contrast table for the primary set. Figure 3 reads this, so it has to be written
# here rather than assembled by hand: keys are "GSE|CellLine|Agent".
write_ranking <- function(nm, out) {
  a <- res[[nm]]
  k <- do.call(rbind, strsplit(a$keys, "|", fixed = TRUE))
  d <- data.frame(GSE = k[,1], Line = k[,2], Agent = k[,3],
                  est = as.numeric(a$obs), null_p = as.numeric(a$null_p),
                  src = a$src, stringsAsFactors = FALSE)
  d$q <- stats::p.adjust(d$null_p, "BH")
  d <- d[order(d$est), ]
  write.csv(d, file.path("results", out), row.names = FALSE)
  cat(sprintf("Wrote results/%s (%d contrasts, %d surviving q<0.10)\n",
              out, nrow(d), sum(d$q < 0.10)))
}
write_ranking("LSC18", "lsc18_singscore_ranking.csv")

cat("\n=== weighted sum versus rank-based enrichment, same contrasts and null ===\n\n")
cat(sprintf("%-8s %5s %10s %9s %11s %9s %9s\n",
            "set","genes","contrasts","beat null","expected","BH q<.10","global P"))
for (nm in names(res)) {
  r <- res[[nm]]
  gp <- (1 + sum(r$global_null >= r$global_obs)) / (length(r$global_null) + 1)
  cat(sprintf("%-8s %5d %10d %9d %11.1f %9d %9.3f\n", nm, r$n, length(r$obs),
              sum(r$null_p < .05), .05*length(r$obs),
              sum(p.adjust(r$null_p,"BH") < .10), gp))
}
cat("\n  (weighted sum, for comparison: 18 of 257 beat null, expected 12.9, 0 survive BH,\n   global P = 0.99 array / 0.39 counts)\n")
