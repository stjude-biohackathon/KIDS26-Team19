#!/usr/bin/env Rscript
# Score every array-based study in the bronze corpus that has a drug and a control arm.
#
# Sequencing studies are excluded here: a GEO series matrix for an RNA-seq submission
# holds metadata only, with counts in per-study supplementary files. They are covered by
# R/score_ncbi_corpus.R, which reads the counts NCBI recomputes for GEO2R-flagged series.
# ARCHS4 or recount3 would extend that further and remain separate work.
#
# Resumable:each study writes its own result file and is skipped on a rerun.
#
#   Rscript R/score_corpus.R [signature]

source(file.path("R", "lsc_scores.R"))
source(file.path("R", "harmonise.R"))

SIG   <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) a[1] else "pLSC6" }
MXDIR <- file.path("bronze", "matrices")
OUT   <- file.path("results", "corpus_scores")
dir.create(MXDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT,   recursive = TRUE, showWarnings = FALSE)

SEQ <- c("GPL24676","GPL18573","GPL16791","GPL11154","GPL20301","GPL9052","GPL10999",
         "GPL15520","GPL21290","GPL17021","GPL19057","GPL13112","GPL11002","GPL18460",
         "GPL25526","GPL26963","GPL30173","GPL29155","GPL20795","GPL21697","GPL24247",
         "GPL19415","GPL11488")

scan <- read.csv(file.path("results","corpus_scan.csv"), stringsAsFactors = FALSE)
man  <- read.csv(file.path("bronze","manifest.csv"),     stringsAsFactors = FALSE)
man  <- man[man$Status == "OK", ]
todo <- scan[scan$usable & !(scan$Platform %in% SEQ), ]
todo$URL <- man$URL[match(todo$GSE, man$SeriesAccession)]
todo <- todo[!is.na(todo$URL), ]
message(sprintf("%d array studies with drug and control arms\n", nrow(todo)))

for (k in seq_len(nrow(todo))) {
  gse <- todo$GSE[k]
  rf  <- file.path(OUT, paste0(gse, ".csv"))
  if (file.exists(rf)) next

  mf <- file.path(MXDIR, paste0(gse, "_series_matrix.txt.gz"))
  if (!file.exists(mf) || file.size(mf) == 0) {
    ok <- system2("curl", c("-sfL","--max-time","600","--retry","3","--retry-delay","4",
                            "-o", shQuote(paste0(mf,".part")), shQuote(todo$URL[k])))
    if (ok != 0) { message(sprintf("[%3d/%d] %-11s download failed", k, nrow(todo), gse)); next }
    file.rename(paste0(mf,".part"), mf)
  }

  res <- try({
    sm <- read_series_matrix(mf)
    if (!nrow(sm$expr) || ncol(sm$expr) < 4) stop("no usable expression table")
    logged <- FALSE
    if (isTRUE(max(sm$expr, na.rm = TRUE) > 50)) { sm$expr <- log2(pmax(sm$expr,0)+1); logged <- TRUE }
    s <- suppressWarnings(score_lsc(sm, SIG))
    h <- harmonise(s$Characteristics)
    d <- cbind(s[, c("GSM","SampleTitle","Score","GSE","Platform","GenesFound","GenesExpected")], h)
    d$Log2Applied <- logged
    d
  }, silent = TRUE)

  if (inherits(res, "try-error")) {
    msg <- substr(sub(".*: ", "", attr(res,"condition")$message), 1, 46)
    message(sprintf("[%3d/%d] %-11s skipped: %s", k, nrow(todo), gse, msg))
    writeLines(msg, file.path(OUT, paste0(gse, ".skip")))
    next
  }
  write.csv(res, rf, row.names = FALSE)
  message(sprintf("[%3d/%d] %-11s scored n=%-3d %d/%d genes", k, nrow(todo), gse,
                  nrow(res), res$GenesFound[1], res$GenesExpected[1]))
}

done <- list.files(OUT, pattern = "\\.csv$", full.names = TRUE)
if (length(done)) {
  all <- do.call(rbind, lapply(done, read.csv, stringsAsFactors = FALSE))
  write.csv(all, file.path("results","corpus_sample_scores.csv"), row.names = FALSE)
  message(sprintf("\n%d studies scored, %d samples -> results/corpus_sample_scores.csv",
                  length(done), nrow(all)))
}
message(sprintf("%d studies skipped", length(list.files(OUT, pattern="\\.skip$"))))
