#!/usr/bin/env Rscript
# Score the RNA-seq series that NCBI publishes counts for.
#
# These are the studies the series-matrix path cannot reach: a sequencing submission's
# series matrix carries metadata only. GEO's GEO2R flag identifies which of them have
# NCBI-recomputed counts, so the target list is knowable in advance rather than by
# downloading and failing.
#
# Resumable: each series writes its own file and is skipped on a rerun.
#
#   Rscript R/score_ncbi_corpus.R [targets.csv]

source(file.path("R", "lsc_scores.R"))
source(file.path("R", "ncbi_counts.R"))
source(file.path("R", "harmonise.R"))

a    <- commandArgs(trailingOnly = TRUE)
TGT  <- if (length(a)) a[1] else file.path("data", "ncbi_counts_targets.csv")
OUT  <- file.path("results", "ncbi_scores")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

sig  <- load_signature("pLSC6")
todo <- read.csv(TGT, stringsAsFactors = FALSE)
message(sprintf("%d series with NCBI counts and a treated/control design\n", nrow(todo)))

for (k in seq_len(nrow(todo))) {
  gse <- todo$GSE[k]
  rf  <- file.path(OUT, paste0(gse, ".csv"))
  if (file.exists(rf) || file.exists(paste0(rf, ".skip"))) next

  res <- try({
    s <- score_ncbi(gse, sig)
    if (all(is.na(s$Characteristics)))
      stop("no sample characteristics (header missing)")
    h <- harmonise(s$Characteristics)
    cbind(s, h)
  }, silent = TRUE)

  if (inherits(res, "try-error")) {
    msg <- substr(sub(".*: ", "", attr(res, "condition")$message), 1, 52)
    message(sprintf("[%2d/%d] %-11s skipped: %s", k, nrow(todo), gse, msg))
    writeLines(msg, paste0(rf, ".skip")); next
  }
  write.csv(res, rf, row.names = FALSE)
  arms <- setdiff(unique(na.omit(res$Agent)), NA)
  message(sprintf("[%2d/%d] %-11s n=%-4d %d/%d genes  arms: %s", k, nrow(todo), gse,
                  nrow(res), res$GenesFound[1], res$GenesExpected[1],
                  if (length(arms)) paste(utils::head(arms, 4), collapse = ", ") else "none"))
}

done <- list.files(OUT, pattern = "\\.csv$", full.names = TRUE)
if (length(done)) {
  all <- do.call(rbind, lapply(done, function(f) {
    d <- read.csv(f, stringsAsFactors = FALSE); d$Genotype <- NULL; d }))
  write.csv(all, file.path("results", "ncbi_sample_scores.csv"), row.names = FALSE)
  message(sprintf("\n%d series scored, %d samples -> results/ncbi_sample_scores.csv",
                  length(done), nrow(all)))
}
message(sprintf("%d skipped", length(list.files(OUT, pattern = "\\.skip$"))))
