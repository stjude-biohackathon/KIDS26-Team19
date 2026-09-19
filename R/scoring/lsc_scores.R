# Leukemic stem cell scoring for GEO series matrix data.
#
# Computes LSC17 (Ng 2016) and pLSC6 (Elsayed 2020) from a GEO series matrix.
# Coefficients live in data/signatures/*.csv, not here, so they can be checked
# against the papers. See data/signatures/README.md for what the scores mean and
# for the standing warning that no coefficient has been verified yet.
#
# Depends only on base R plus curl on the PATH. No Bioconductor.

SIG_DIR         <- file.path("data", "signatures")
ANNOT_CACHE_DIR <- file.path("bronze", "platform-annot")

SIGNATURES <- list(
  LSC17 = list(file = "lsc17_ng2016.csv",
               cite = "Ng et al. Nature 2016;540:433-437"),
  pLSC6 = list(file = "plsc6_elsayed2020.csv",
               cite = "Elsayed et al. Leukemia 2020;34:735-745")
)

# ---- signatures ------------------------------------------------------------

load_signature <- function(name, sig_dir = SIG_DIR) {
  if (!name %in% names(SIGNATURES)) {
    stop("Unknown signature '", name, "'. Available: ",
         paste(names(SIGNATURES), collapse = ", "), call. = FALSE)
  }
  path <- file.path(sig_dir, SIGNATURES[[name]]$file)
  if (!file.exists(path)) stop("Signature file not found: ", path, call. = FALSE)

  sig <- read.csv(path, stringsAsFactors = FALSE, na.strings = "")
  stopifnot(all(c("gene", "coefficient", "aliases", "verified") %in% names(sig)))
  if (anyNA(sig$coefficient)) stop("Missing coefficient in ", path, call. = FALSE)

  sig$aliases  <- ifelse(is.na(sig$aliases), "", sig$aliases)
  sig$verified <- toupper(as.character(sig$verified)) %in% c("TRUE", "YES", "1")
  attr(sig, "name") <- name
  attr(sig, "cite") <- SIGNATURES[[name]]$cite
  sig
}

# Every symbol that should be treated as this gene, primary name first.
symbols_for <- function(row) {
  alt <- unlist(strsplit(row$aliases, "[;,|]"))
  unique(trimws(c(row$gene, alt[nzchar(trimws(alt))])))
}

# ---- platform annotation ---------------------------------------------------

# GEO publishes a curated annotation table for some platforms and not others.
# Two sources, in order:
#   1. <GPL>.annot.gz, which has a proper "Gene symbol" column.
#   2. the platform table from acc.cgi, where the symbol column is named
#      differently on every array family (ORF, GENE_SYMBOL, GENE, ...).
# Both are cached on first use.
read_geo_table <- function(con) {
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("No platform table found.", call. = FALSE)
    if (grepl("^!platform_table_begin", line)) break
  }
  tab <- read.delim(con, quote = "", comment.char = "", stringsAsFactors = FALSE,
                    check.names = FALSE, na.strings = c("", "NA"), fill = TRUE)
  tab[!grepl("^!platform_table_end", tab[[1]]), , drop = FALSE]
}

# Pick the column that actually contains gene symbols. Column names are not
# reliable across array families, so choose by overlap with the genes we need and
# fall back to name priority when nothing overlaps.
pick_symbol_column <- function(tab, prefer_symbols = NULL) {
  nm   <- names(tab)
  keep <- seq_along(nm)[-1]
  keep <- keep[!grepl("^(SEQUENCE|SPOT_ID|GB_ACC|RANGE_|PROBE)", nm[keep], ignore.case = TRUE)]
  if (!length(keep)) keep <- seq_along(nm)[-1]

  if (length(prefer_symbols)) {
    hits <- vapply(keep, function(j) {
      v <- as.character(tab[[j]])
      length(intersect(unique(v[!is.na(v)]), prefer_symbols))
    }, integer(1))
    if (length(hits) && max(hits) > 0) return(keep[which.max(hits)])
  }
  for (pat in c("^Gene symbol$", "^GENE_SYMBOL$", "symbol", "^ORF$", "^GENE$")) {
    j <- grep(pat, nm, ignore.case = TRUE)
    j <- j[j != 1]
    if (length(j)) return(j[1])
  }
  NA_integer_
}

fetch_platform_annotation <- function(gpl, cache_dir = ANNOT_CACHE_DIR,
                                      prefer_symbols = NULL) {
  if (!grepl("^GPL\\d+$", gpl)) stop("Not a platform accession: ", gpl, call. = FALSE)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  gz   <- file.path(cache_dir, paste0(gpl, ".annot.gz"))
  soft <- file.path(cache_dir, paste0(gpl, ".table.txt"))
  stub <- paste0(substr(gpl, 1, nchar(gpl) - 3), "nnn")

  get <- function(url, dest) {
    part <- paste0(dest, ".part")
    ok <- system2("curl", c("-sfL", "--max-time", "300", "-o", shQuote(part), shQuote(url)))
    if (ok != 0 || !file.exists(part) || file.size(part) == 0) { unlink(part); return(FALSE) }
    file.rename(part, dest); TRUE
  }

  tab <- NULL
  if (file.exists(gz) && file.size(gz) > 0) {
    con <- gzfile(gz, "rt"); tab <- try(read_geo_table(con), silent = TRUE); close(con)
    if (inherits(tab, "try-error")) tab <- NULL
  }
  if (is.null(tab)) {
    url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/platforms/%s/%s/annot/%s.annot.gz",
                   stub, gpl, gpl)
    if (get(url, gz)) {
      con <- gzfile(gz, "rt"); tab <- try(read_geo_table(con), silent = TRUE); close(con)
      if (inherits(tab, "try-error")) tab <- NULL
    }
  }
  if (is.null(tab)) {
    message("  no annot file for ", gpl, "; falling back to the platform table")
    if (!file.exists(soft) || file.size(soft) == 0) {
      url <- sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&targ=self&view=data&form=text", gpl)
      if (!get(url, soft)) stop("Could not retrieve any annotation for ", gpl,
                                ". Supply probe_map= instead.", call. = FALSE)
    }
    con <- file(soft, "rt"); tab <- read_geo_table(con); close(con)
  }

  sym_col <- pick_symbol_column(tab, prefer_symbols)
  if (!is.na(sym_col)) {
    return(data.frame(probe  = as.character(tab[[1]]),
                      symbol = as.character(tab[[sym_col]]),
                      stringsAsFactors = FALSE))
  }

  # Affymetrix ST, Clariom and HTA arrays carry no symbol column. The symbol sits
  # inside an assignment string, as "NM_198317 // RefSeq // ... (KLHL17), mRNA.".
  # Selected by position: these tables often repeat a column name such as SPOT_ID,
  # and lookup by name silently returns the wrong one. Heuristic, so the coverage
  # check in score_lsc() is the real guard.
  is_assign <- vapply(seq_along(tab), function(j)
    isTRUE(mean(grepl(" // ", as.character(tab[[j]]), fixed = TRUE), na.rm = TRUE) > 0.3),
    logical(1))
  if (any(is_assign)) {
    message("  no symbol column for ", gpl, "; parsing Affymetrix assignments")
    v   <- as.character(tab[[which(is_assign)[1]]])
    pat <- "\\(([A-Za-z0-9][A-Za-z0-9._-]{0,14})\\)[,.]?[[:space:]]*(mRNA|non-coding|transcript)"
    m   <- regexpr(pat, v)
    sym <- rep(NA_character_, length(v))
    sym[m > 0] <- sub(paste0("^", pat, "$"), "\\1", regmatches(v, m))
    if (isTRUE(mean(!is.na(sym)) > 0.05)) {
      return(data.frame(probe = as.character(tab[[1]]), symbol = sym,
                        stringsAsFactors = FALSE))
    }
  }

  stop("No gene symbol column in the annotation for ", gpl,
       ". Supply probe_map= instead.", call. = FALSE)
}

# ---- series matrix ---------------------------------------------------------

# Returns the expression matrix plus the !Sample_* metadata lines, which carry
# the cell line and treatment fields the project keys on.
read_series_matrix <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)

  meta <- list()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("No expression table in ", path, call. = FALSE)
    if (grepl("^!series_matrix_table_begin", line)) break
    if (grepl("^!(Sample|Series)_", line)) {
      parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
      key   <- sub("^!", "", parts[1])
      val   <- gsub('^"|"$', "", parts[-1])
      meta[[key]] <- c(meta[[key]], list(val))
    }
  }

  tab <- read.delim(con, quote = "\"", comment.char = "", stringsAsFactors = FALSE,
                    check.names = FALSE, na.strings = c("", "NA", "null"))
  tab <- tab[!grepl("^!series_matrix_table_end", tab[[1]]), , drop = FALSE]

  probes <- as.character(tab[[1]])
  expr   <- as.matrix(tab[, -1, drop = FALSE])
  storage.mode(expr) <- "double"
  rownames(expr) <- probes

  gpl <- if (!is.null(meta$Series_platform_id)) meta$Series_platform_id[[1]][1] else NA_character_
  gse <- if (!is.null(meta$Series_geo_accession)) meta$Series_geo_accession[[1]][1] else NA_character_

  list(expr = expr, meta = meta, gpl = gpl, gse = gse)
}

# Flatten repeated !Sample_characteristics_ch1 lines into one row per sample.
sample_table <- function(sm) {
  gsm <- sm$meta$Sample_geo_accession[[1]]
  out <- data.frame(GSM = gsm, stringsAsFactors = FALSE)
  if (!is.null(sm$meta$Sample_title)) out$SampleTitle <- sm$meta$Sample_title[[1]]
  chars <- sm$meta$Sample_characteristics_ch1
  if (!is.null(chars)) {
    out$Characteristics <- apply(
      do.call(rbind, lapply(chars, function(v) {
        length(v) <- length(gsm); ifelse(is.na(v), "", v)
      })), 2, function(col) paste(col[nzchar(col)], collapse = " | ")
    )
  }
  out
}

# ---- probe to gene ---------------------------------------------------------

# Collapse probes to genes. Several probes map to one gene on most arrays;
# taking the probe with the highest mean expression is the usual convention and
# is recorded in the result so the choice stays visible.
collapse_to_genes <- function(expr, probe_map, wanted) {
  pm <- probe_map[!is.na(probe_map$probe) & !is.na(probe_map$symbol), , drop = FALSE]
  pm <- pm[pm$probe %in% rownames(expr) & nzchar(pm$symbol), , drop = FALSE]
  # GEO packs multi-gene probes as "A///B"; such probes are ambiguous, so drop them.
  pm <- pm[!grepl("///", pm$symbol, fixed = TRUE), , drop = FALSE]

  hits <- pm[pm$symbol %in% wanted, , drop = FALSE]
  if (!nrow(hits)) return(list(mat = NULL, probes = hits[0, ], method = "max_mean"))

  sub  <- expr[hits$probe, , drop = FALSE]
  mean_expr <- rowMeans(sub, na.rm = TRUE)
  best <- tapply(seq_len(nrow(hits)), hits$symbol, function(i) i[which.max(mean_expr[i])])
  best <- unlist(best)

  mat <- sub[best, , drop = FALSE]
  rownames(mat) <- hits$symbol[best]
  list(mat = mat, probes = hits[best, , drop = FALSE], method = "max_mean")
}

# ---- scoring ---------------------------------------------------------------

#' Score samples with an LSC signature.
#'
#' @param sm            result of read_series_matrix()
#' @param signature     "LSC17" or "pLSC6"
#' @param probe_map     optional data.frame(probe, symbol); fetched from GEO if NULL
#' @param min_coverage  minimum fraction of signature genes that must be found
#' @param allow_unverified  proceed although coefficients are unchecked
#' @param sig_dir       directory holding the signature CSVs
score_lsc <- function(sm, signature = "LSC17", probe_map = NULL,
                      min_coverage = 0.8, allow_unverified = FALSE,
                      sig_dir = SIG_DIR) {
  sig <- load_signature(signature, sig_dir = sig_dir)

  unverified <- sum(!sig$verified)
  if (unverified > 0) {
    msg <- sprintf(
      "%s: %d of %d coefficients are unverified against %s. Scores are a plumbing test, not a result.",
      signature, unverified, nrow(sig), attr(sig, "cite"))
    if (!allow_unverified) {
      stop(msg, "\n  Verify ", file.path(sig_dir, SIGNATURES[[signature]]$file),
           " and set verified=TRUE, or pass allow_unverified=TRUE.", call. = FALSE)
    }
    warning(msg, call. = FALSE, immediate. = TRUE)
  }

  wanted <- unique(unlist(lapply(seq_len(nrow(sig)), function(i) symbols_for(sig[i, ]))))

  if (is.null(probe_map)) {
    if (is.na(sm$gpl)) stop("No platform id in the series matrix; supply probe_map=.", call. = FALSE)
    probe_map <- fetch_platform_annotation(sm$gpl, prefer_symbols = wanted)
  }

  coll    <- collapse_to_genes(sm$expr, probe_map, wanted)
  if (is.null(coll$mat)) stop("None of the signature genes were found on ", sm$gpl, call. = FALSE)

  # Map whichever alias the platform used back onto the signature's primary gene.
  found <- rownames(coll$mat)
  owner <- vapply(found, function(s) {
    hit <- which(vapply(seq_len(nrow(sig)), function(i) s %in% symbols_for(sig[i, ]), logical(1)))
    if (length(hit)) sig$gene[hit[1]] else NA_character_
  }, character(1))

  mat <- coll$mat[!is.na(owner), , drop = FALSE]
  rownames(mat) <- owner[!is.na(owner)]
  mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]

  coverage <- nrow(mat) / nrow(sig)
  missing  <- setdiff(sig$gene, rownames(mat))
  if (coverage < min_coverage) {
    stop(sprintf("Only %d/%d %s genes found on %s (%.0f%%, need %.0f%%). Missing: %s",
                 nrow(mat), nrow(sig), signature, sm$gpl, 100 * coverage,
                 100 * min_coverage, paste(missing, collapse = ", ")), call. = FALSE)
  }

  # The score is a weighted sum, so an un-logged matrix produces values on a wildly
  # different scale. Both published scores were derived on log-scale expression.
  scale_note <- "log-like"
  mx <- suppressWarnings(max(mat, na.rm = TRUE))
  if (is.finite(mx) && mx > 50) {
    scale_note <- "NOT log-scale"
    warning(sprintf(
      "%s: max signature-gene value is %.1f, so this matrix looks un-logged. Both scores assume log-scale input; consider log2(x + 1) first.",
      sm$gse, mx), call. = FALSE, immediate. = TRUE)
  }

  coefs <- setNames(sig$coefficient, sig$gene)[rownames(mat)]
  raw   <- as.numeric(crossprod(mat, coefs))

  # Raw scores carry the units of the input matrix, so they are meaningless across
  # studies. The within-study z is the comparable quantity.
  z <- if (stats::sd(raw, na.rm = TRUE) > 0) as.numeric(scale(raw)) else rep(NA_real_, length(raw))

  out <- sample_table(sm)
  out$Signature     <- signature
  out$Score         <- raw
  out$ScoreZ        <- z
  out$GenesFound    <- nrow(mat)
  out$GenesExpected <- nrow(sig)
  out$Coverage      <- round(coverage, 3)
  out$MissingGenes  <- paste(missing, collapse = ";")
  out$GSE           <- sm$gse
  out$Platform      <- sm$gpl
  out$CollapseMethod    <- coll$method
  out$InputScale        <- scale_note
  out$CoefficientsVerified <- unverified == 0

  attr(out, "genes")  <- rownames(mat)
  attr(out, "probes") <- coll$probes
  out
}
