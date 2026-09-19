#!/usr/bin/env Rscript
# Scoring from NCBI-generated RNA-seq counts.
#
# A GEO series matrix for a sequencing submission holds metadata only, which puts most
# recent studies out of reach of the series-matrix path in R/lsc_scores.R. For human and
# mouse RNA-seq, however, NCBI recomputes counts uniformly and publishes them per series.
# GEO's own GEO2R flag says when they exist, so availability is knowable before download.
#
# Rows of those files are Entrez GeneIDs, not symbols, so the signature tables carry a
# verified `entrez` column and this module keys on it.
#
# Normalisation: counts are converted to counts per million and log2(CPM + 1) is taken.
# The published score was fitted on log2(RPKM + 1), which additionally divides by gene
# length. That difference is a constant per gene, and every contrast here is computed
# within a single study where all samples share it, so it cancels in the treated-minus-
# control difference and in the within-study z-score. It would NOT cancel in a comparison
# of absolute scores across studies, which this pipeline never makes.

COUNTS_CACHE <- file.path("bronze", "ncbi-counts")
GEO_DL <- "https://www.ncbi.nlm.nih.gov/geo/download/"

# ---- retrieval -------------------------------------------------------------

# NCBI publishes counts against a fixed assembly build; the file name encodes it.
ncbi_counts_url <- function(gse, build = "GRCh38.p13") {
  sprintf("%s?type=rnaseq_counts&acc=%s&format=file&file=%s_raw_counts_%s_NCBI.tsv.gz",
          GEO_DL, gse, gse, build)
}

fetch_ncbi_counts <- function(gse, cache_dir = COUNTS_CACHE, build = "GRCh38.p13") {
  if (!grepl("^GSE\\d+$", gse)) stop("Not a series accession: ", gse, call. = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  fp <- file.path(cache_dir, sprintf("%s_raw_counts_%s.tsv.gz", gse, build))

  if (!file.exists(fp) || file.size(fp) == 0) {
    part <- paste0(fp, ".part")
    ok <- system2("curl", c("-sfL", "--max-time", "600", "--retry", "2",
                            "-o", shQuote(part), shQuote(ncbi_counts_url(gse, build))))
    if (ok != 0 || !file.exists(part) || file.size(part) < 1000) {
      unlink(part)
      stop("No NCBI-generated counts for ", gse,
           " (GEO2R must report them; check ExpressionDataAvailability).", call. = FALSE)
    }
    # A 404 is served as an HTML page, which is not gzip.
    if (!identical(as.integer(readBin(part, "raw", 2)), c(31L, 139L))) {
      unlink(part); stop("Response for ", gse, " is not a gzip file.", call. = FALSE)
    }
    file.rename(part, fp)
  }
  fp
}

# ---- matrix ----------------------------------------------------------------

#' Signature-gene expression for one series, from NCBI counts.
#'
#' Returns a list shaped like read_series_matrix() so the rest of the pipeline can
#' consume it unchanged: `expr` carries log2(CPM + 1) with gene SYMBOLS as row names.
#'
#' @param gse   series accession
#' @param sig   signature table from load_signature(); must carry an `entrez` column
#' @param keep_all  return every gene rather than only the signature (used by the
#'                  permutation tests, which need a sampling pool)
ncbi_counts_matrix <- function(gse, sig, cache_dir = COUNTS_CACHE, keep_all = FALSE) {
  if (!"entrez" %in% names(sig) || all(is.na(sig$entrez) | sig$entrez == ""))
    stop("Signature has no Entrez ids; the counts path needs them.", call. = FALSE)

  fp  <- fetch_ncbi_counts(gse, cache_dir)
  con <- gzfile(fp, "rt"); on.exit(close(con), add = TRUE)
  tab <- utils::read.delim(con, check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(tab) || ncol(tab) < 3) stop("Counts table for ", gse, " is unusable.", call. = FALSE)

  ids <- as.character(tab[[1]])
  m   <- as.matrix(tab[, -1, drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- ids

  # Library-size normalisation must use the WHOLE library, not the signature subset.
  lib <- colSums(m, na.rm = TRUE)
  if (any(lib <= 0)) stop("Zero library size in ", gse, call. = FALSE)
  cpm <- log2(sweep(m, 2, lib, "/") * 1e6 + 1)

  if (keep_all) {
    out <- cpm[stats::complete.cases(cpm), , drop = FALSE]
  } else {
    want <- as.character(sig$entrez)
    hit  <- match(want, rownames(cpm))
    if (all(is.na(hit))) stop("No signature genes found in counts for ", gse, call. = FALSE)
    out <- cpm[hit[!is.na(hit)], , drop = FALSE]
    rownames(out) <- sig$gene[!is.na(hit)]
  }

  list(expr = out, meta = list(), gpl = "NCBI-counts", gse = gse,
       lib_size = lib, source = "ncbi_rnaseq_counts")
}

# ---- sample metadata -------------------------------------------------------

# The counts file carries no metadata, so characteristics still come from the series
# matrix header, which exists for sequencing submissions even when the table does not.
ncbi_sample_metadata <- function(gse, header_dir = file.path("bronze", "headers")) {
  fp <- file.path(header_dir, paste0(gse, ".hdr"))
  if (!file.exists(fp)) return(NULL)
  ln  <- readLines(fp, warn = FALSE)
  get <- function(key) {
    i <- grep(paste0("^!", key, "\t"), ln)
    if (!length(i)) return(NULL)
    lapply(i, function(j) gsub('^"|"$', "", strsplit(ln[j], "\t", fixed = TRUE)[[1]][-1]))
  }
  gsm <- get("Sample_geo_accession"); if (is.null(gsm)) return(NULL)
  gsm <- gsm[[1]]
  ch  <- get("Sample_characteristics_ch1")
  chars <- if (is.null(ch)) rep("", length(gsm)) else
    apply(do.call(rbind, lapply(ch, function(v) { length(v) <- length(gsm)
                                                  ifelse(is.na(v), "", v) })),
          2, function(col) paste(col[nzchar(col)], collapse = " | "))
  ttl <- get("Sample_title")
  data.frame(GSM = gsm,
             SampleTitle = if (is.null(ttl)) NA_character_ else ttl[[1]],
             Characteristics = chars, stringsAsFactors = FALSE)
}

# ---- scoring ---------------------------------------------------------------

#' Score one series from NCBI counts.
score_ncbi <- function(gse, sig, cache_dir = COUNTS_CACHE, min_coverage = 1.0,
                       allow_unverified = FALSE) {
  if (any(!sig$verified) && !allow_unverified)
    stop("Signature coefficients are unverified.", call. = FALSE)

  sm  <- ncbi_counts_matrix(gse, sig, cache_dir)
  cov <- nrow(sm$expr) / nrow(sig)
  if (cov < min_coverage)
    stop(sprintf("Only %d/%d signature genes in counts for %s (%.0f%%).",
                 nrow(sm$expr), nrow(sig), gse, 100 * cov), call. = FALSE)

  co    <- setNames(sig$coefficient, sig$gene)[rownames(sm$expr)]
  score <- as.numeric(t(sm$expr) %*% co)

  out <- data.frame(GSM = colnames(sm$expr), Score = score, GSE = gse,
                    Platform = "NCBI-counts", GenesFound = nrow(sm$expr),
                    GenesExpected = nrow(sig), stringsAsFactors = FALSE)
  md <- ncbi_sample_metadata(gse)
  if (!is.null(md)) out <- merge(out, md, by = "GSM", all.x = TRUE, sort = FALSE)
  else out$Characteristics <- NA_character_
  out$LibSize <- sm$lib_size[match(out$GSM, names(sm$lib_size))]
  out
}
