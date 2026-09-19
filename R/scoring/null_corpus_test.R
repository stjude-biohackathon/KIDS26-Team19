#!/usr/bin/env Rscript
# Permutation test across the corpus: do the 95 contrasts depend on pLSC6 specifically,
# or would any six genes with these weights produce them?
#
# Run on the 37 corpus studies (supersedes the 4-study pilot test, see git history). Each study's
# null draws the same number of genes the observed signature had in that study, so
# 5-gene studies are compared against 5-gene nulls.
#
#   Rscript R/null_corpus_test.R [nperm]

set.seed(20260918)
source(file.path("R", "lsc_scores.R"))
NPERM <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) as.integer(a[1]) else 1000L }

sig    <- load_signature("pLSC6"); genes <- sig$gene; coefs <- sig$coefficient
wanted <- unique(unlist(lapply(seq_len(nrow(sig)), function(i) symbols_for(sig[i, ]))))

ss <- read.csv(file.path("results", "corpus_sample_scores.csv"), stringsAsFactors = FALSE)
ss$Line <- ifelse(!is.na(ss$CellLine), ss$CellLine,
            ifelse(!is.na(ss$Genotype), ss$Genotype, ss$GSE))

studies <- list()
for (gse in unique(ss$GSE)) {
  mf <- file.path("bronze", "matrices", paste0(gse, "_series_matrix.txt.gz"))
  if (!file.exists(mf)) next
  sm <- try(read_series_matrix(mf), silent = TRUE); if (inherits(sm, "try-error")) next
  if (isTRUE(max(sm$expr, na.rm = TRUE) > 50)) sm$expr <- log2(pmax(sm$expr, 0) + 1)
  pm <- try(fetch_platform_annotation(sm$gpl, prefer_symbols = wanted), silent = TRUE)
  if (inherits(pm, "try-error")) next

  pmk <- pm[!is.na(pm$probe) & !is.na(pm$symbol), ]
  pmk <- pmk[pmk$probe %in% rownames(sm$expr) & nzchar(pmk$symbol) &
             !grepl("///", pmk$symbol, fixed = TRUE), ]
  if (nrow(pmk) < 500) next
  sub <- sm$expr[pmk$probe, , drop = FALSE]
  mu  <- rowMeans(sub, na.rm = TRUE)
  best <- unlist(tapply(seq_len(nrow(pmk)), pmk$symbol, function(i) i[which.max(mu[i])]))
  M <- sub[best, , drop = FALSE]; rownames(M) <- pmk$symbol[best]
  M <- M[stats::complete.cases(M), , drop = FALSE]
  for (i in seq_len(nrow(sig))) {
    syn <- symbols_for(sig[i, ]); hit <- which(rownames(M) %in% syn)
    if (length(hit) && !(sig$gene[i] %in% rownames(M))) rownames(M)[hit[1]] <- sig$gene[i]
  }
  have <- genes[genes %in% rownames(M)]; if (length(have) < 4) next

  d <- ss[ss$GSE == gse, ]; d <- d[match(colnames(M), d$GSM), ]
  grp <- list()
  for (ln in unique(na.omit(d$Line))) {
    i  <- !is.na(d$Line) & d$Line == ln
    ct <- which(i & !is.na(d$Agent) & d$Agent == "control"); if (length(ct) < 2) next
    for (a in setdiff(unique(na.omit(d$Agent[i])), "control")) {
      tr <- which(i & !is.na(d$Agent) & d$Agent == a); if (length(tr) < 2) next
      grp[[length(grp) + 1]] <- list(key = paste(gse, ln, a, sep = "|"),
                                     agent = a, line = ln, tr = tr, ct = ct)
    }
  }
  if (!length(grp)) next
  studies[[gse]] <- list(M = M, grp = grp, have = have,
                         co = setNames(coefs, genes)[have],
                         pool = which(apply(M, 1, stats::sd) > 0))
}
message(sprintf("%d studies loaded, %d contrasts", length(studies),
                sum(vapply(studies, function(s) length(s$grp), integer(1)))))

eff <- function(M, idx, grp, co) {
  z <- as.numeric(scale(as.numeric(t(M[idx, , drop = FALSE]) %*% co)))
  vapply(grp, function(g) mean(z[g$tr]) - mean(z[g$ct]), numeric(1))
}
obs  <- unlist(lapply(studies, function(s) eff(s$M, match(s$have, rownames(s$M)), s$grp, s$co)))
keys <- unlist(lapply(studies, function(s) vapply(s$grp, function(g) g$key, character(1))))

message(sprintf("Permuting %d random signatures ...", NPERM))
null <- matrix(NA_real_, NPERM, length(obs), dimnames = list(NULL, keys))
for (p in seq_len(NPERM))
  null[p, ] <- unlist(lapply(studies, function(s)
                  eff(s$M, sample(s$pool, length(s$have)), s$grp, s$co)))

emp <- vapply(seq_along(obs), function(j)
  (1 + sum(abs(null[, j]) >= abs(obs[j]))) / (NPERM + 1), numeric(1))
parts <- do.call(rbind, strsplit(keys, "|", fixed = TRUE))
res <- data.frame(GSE = parts[,1], Line = parts[,2], Agent = parts[,3],
                  observed = round(obs, 3), null_p = round(emp, 4),
                  stringsAsFactors = FALSE)
res$null_q <- round(p.adjust(res$null_p, "BH"), 4)
res <- res[order(res$null_p), ]
# Pinometostat genotype split, the one cell-line result that survived at pilot scale
MLLR <- c("NOMO-1","THP-1","MV4-11","MOLM-13","MLL rearranged")
MLLW <- c("HL-60","U-937","OCI-AML3","MLL wild type")
pin  <- grep("[|]Pinometostat$", keys)
pl   <- sub("^[^|]+[|]", "", sub("[|]Pinometostat$", "", keys[pin]))
jr <- pin[pl %in% MLLR]; jw <- pin[pl %in% MLLW]
sep_obs <- if (length(jr) && length(jw)) mean(obs[jr]) - mean(obs[jw]) else NA_real_
sep_null <- if (length(jr) && length(jw))
  rowMeans(null[, jr, drop=FALSE]) - rowMeans(null[, jw, drop=FALSE]) else NULL
if (!is.na(sep_obs))
  cat(sprintf("\nPinometostat, MLL-rearranged minus wild type: %.3f  (P = %.4f)\n",
              sep_obs, (1 + sum(sep_null <= sep_obs)) / (NPERM + 1)))

saveRDS(list(obs = obs, null = null, keys = keys,
             global_obs = mean(abs(obs)), global_null = rowMeans(abs(null)),
             sep_obs = sep_obs, sep_null = sep_null),
        file.path("results", "corpus_null_draws.rds"))
write.csv(res, file.path("results", "corpus_null_test.csv"), row.names = FALSE)

cat(sprintf("\n=== %d contrasts vs %d random signatures ===\n\n", nrow(res), NPERM))
cat(sprintf("  beating the null, P<0.05 : %d (%.0f%%)\n", sum(res$null_p<.05), 100*mean(res$null_p<.05)))
cat(sprintf("  after BH, q<0.10         : %d\n", sum(res$null_q<.10)))
cat(sprintf("  global mean |effect|     : observed %.3f vs random %.3f (p = %.4f)\n\n",
            mean(abs(obs)), mean(abs(null)),
            (1 + sum(rowMeans(abs(null)) >= mean(abs(obs)))) / (NPERM + 1)))
cat("  most signature-specific contrasts:\n\n")
print(head(data.frame(Agent = substr(res$Agent,1,26), Line = substr(res$Line,1,18),
                      est = res$observed, null_p = res$null_p, null_q = res$null_q), 18),
      row.names = FALSE)
cat("\nWrote results/corpus_null_test.csv\n")
