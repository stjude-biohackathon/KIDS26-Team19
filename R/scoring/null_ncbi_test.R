#!/usr/bin/env Rscript
# Permutation test for the NCBI-counts contrasts.
#
# Same question as R/null_corpus_test.R: do these effects depend on the pLSC6 genes, or
# would any six genes with these weights produce them? Processes one study at a time and
# discards each expression matrix after use, because the counts matrices are large
# (up to ~40,000 genes x 300 samples).
#
#   Rscript R/null_ncbi_test.R [nperm]

set.seed(20260918)
source(file.path("R", "lsc_scores.R"))
source(file.path("R", "ncbi_counts.R"))
NPERM <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) as.integer(a[1]) else 1000L }

sig <- load_signature("pLSC6")
ss  <- read.csv(file.path("results", "ncbi_sample_scores.csv"), stringsAsFactors = FALSE)
ss$Line <- ifelse(!is.na(ss$CellLine), ss$CellLine, ss$GSE)

obs_all <- numeric(0); keys_all <- character(0)
null_cols <- list()

for (gse in unique(ss$GSE)) {
  d <- ss[ss$GSE == gse, ]
  if (all(is.na(d$Agent))) next

  M <- try(ncbi_counts_matrix(gse, sig, keep_all = TRUE)$expr, silent = TRUE)
  if (inherits(M, "try-error")) { message("  ", gse, " matrix unavailable"); next }
  d <- d[match(colnames(M), d$GSM), ]

  grp <- list()
  for (ln in unique(na.omit(d$Line))) {
    i  <- !is.na(d$Line) & d$Line == ln
    ct <- which(i & !is.na(d$Agent) & d$Agent == "control"); if (length(ct) < 2) next
    for (a in setdiff(unique(na.omit(d$Agent[i])), "control")) {
      tr <- which(i & !is.na(d$Agent) & d$Agent == a); if (length(tr) < 2) next
      grp[[length(grp) + 1]] <- list(key = paste(gse, ln, a, sep = "|"), tr = tr, ct = ct)
    }
  }
  if (!length(grp)) { rm(M); gc(verbose = FALSE); next }

  idx  <- match(as.character(sig$entrez), rownames(M))
  # rownames are Entrez ids here; the signature rows were renamed on the scoring path
  if (anyNA(idx)) idx <- match(sig$gene, rownames(M))
  if (anyNA(idx)) { message("  ", gse, " signature genes absent"); rm(M); gc(verbose=FALSE); next }

  pool <- which(apply(M, 1, stats::sd) > 0)
  eff  <- function(rows) {
    z <- as.numeric(scale(as.numeric(t(M[rows, , drop = FALSE]) %*% sig$coefficient)))
    vapply(grp, function(g) mean(z[g$tr]) - mean(z[g$ct]), numeric(1))
  }
  o <- eff(idx)
  nl <- vapply(seq_len(NPERM), function(p) eff(sample(pool, nrow(sig))),
               numeric(length(grp)))
  if (length(grp) == 1L) nl <- matrix(nl, nrow = 1)

  obs_all  <- c(obs_all, o)
  keys_all <- c(keys_all, vapply(grp, function(g) g$key, character(1)))
  null_cols[[length(null_cols) + 1]] <- t(nl)
  message(sprintf("  %-11s %d genes, %d contrasts", gse, nrow(M), length(grp)))
  rm(M); gc(verbose = FALSE)
}

null <- do.call(cbind, null_cols)
colnames(null) <- keys_all
emp <- vapply(seq_along(obs_all), function(j)
  (1 + sum(abs(null[, j]) >= abs(obs_all[j]))) / (NPERM + 1), numeric(1))

parts <- do.call(rbind, strsplit(keys_all, "|", fixed = TRUE))
res <- data.frame(GSE = parts[,1], Line = parts[,2], Agent = parts[,3],
                  observed = round(obs_all, 3), null_p = round(emp, 4),
                  stringsAsFactors = FALSE)
res$null_q <- round(p.adjust(res$null_p, "BH"), 4)
res <- res[order(res$null_p), ]
write.csv(res, file.path("results", "ncbi_null_test.csv"), row.names = FALSE)
saveRDS(list(obs = obs_all, null = null, keys = keys_all,
             global_obs = mean(abs(obs_all)), global_null = rowMeans(abs(null))),
        file.path("results", "ncbi_null_draws.rds"))

cat(sprintf("\n=== %d contrasts vs %d random signatures ===\n\n", nrow(res), NPERM))
cat(sprintf("  beat the null, P<0.05 : %d (%.1f%%), chance expectation %.1f\n",
            sum(res$null_p < .05), 100*mean(res$null_p < .05), .05*nrow(res)))
cat(sprintf("  after BH, q<0.10      : %d\n", sum(res$null_q < .10)))
cat(sprintf("  global mean |effect|  : observed %.3f vs random %.3f (P = %.4f)\n\n",
            mean(abs(obs_all)), mean(abs(null)),
            (1 + sum(rowMeans(abs(null)) >= mean(abs(obs_all)))) / (NPERM + 1)))
print(utils::head(data.frame(Agent = substr(res$Agent,1,24), Line = substr(res$Line,1,16),
                             est = res$observed, null_p = res$null_p, null_q = res$null_q), 14),
      row.names = FALSE)
