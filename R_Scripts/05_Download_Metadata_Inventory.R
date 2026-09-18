# Standalone GEO discovery, matrix download, platform resolution and CSV reports.
# Rscript R_Scripts/05_Download_Metadata_Inventory.R
# Rscript R_Scripts/05_Download_Metadata_Inventory.R input.xlsx output_directory
# Packages: xml2; readxl for Excel; GEOquery only for missing-platform fallback.
# Sourcing this file defines functions without starting downloads.
# NCBI workflow: https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html#FTP

.geo_empty <- function() {
  data.frame(SeriesAccession = character(), Platform = character(),
             MatrixFile = character(), URL = character(), LocalFile = character(),
             Status = character(), ErrorMessage = character(),
             stringsAsFactors = FALSE)
}

.geo_row <- function(gse, status, error = NA_character_, file = NA_character_,
                     url = NA_character_, platform = NA_character_) {
  data.frame(SeriesAccession = gse, Platform = platform, MatrixFile = file,
             URL = url, LocalFile = NA_character_, Status = status,
             ErrorMessage = error, stringsAsFactors = FALSE)
}

.geo_controls <- function(retries, delay_seconds, timeout_seconds) {
  stopifnot(length(retries) == 1L, is.finite(retries), retries >= 0,
            retries == floor(retries), length(delay_seconds) == 1L,
            is.finite(delay_seconds), delay_seconds >= 0,
            length(timeout_seconds) == 1L, is.finite(timeout_seconds),
            timeout_seconds > 0)
}

# retries counts additional attempts after the first request.
.geo_retry <- function(action, retries, delay_seconds) {
  for (attempt in seq_len(retries + 1L)) {
    result <- tryCatch(list(value = action()), error = function(e) list(error = e))
    if (is.null(result$error)) return(result$value)
    if (attempt > retries) stop(result$error)
    Sys.sleep(delay_seconds * 2^(attempt - 1L))
  }
}

.geo_fetch <- function(url, destination) {
  status <- suppressWarnings(utils::download.file(
    url, destination, mode = "wb", quiet = TRUE, method = "libcurl"
  ))
  if (!identical(status, 0L)) stop("Download failed: ", url)
  invisible(destination)
}

.geo_listing <- function(url) {
  file <- tempfile(fileext = ".html")
  on.exit(unlink(file))
  .geo_fetch(url, file)
  xml2::xml_attr(xml2::xml_find_all(xml2::read_html(file), "//a[@href]"), "href")
}

# Read the actual directory listing: handles both GSE..._series_matrix.txt.gz
# and GSE...-GPL..._series_matrix.txt.gz, without guessing platform IDs.
discover_geo_matrices <- function(gse_ids, retries = 2L, delay_seconds = 0.4,
                                  timeout_seconds = 120) {
  if (!requireNamespace("xml2", quietly = TRUE)) stop("Install the xml2 package first.")
  .geo_controls(retries, delay_seconds, timeout_seconds)
  old <- options(timeout = max(getOption("timeout", 60), timeout_seconds))
  on.exit(options(old))
  ids <- unique(toupper(trimws(as.character(gse_ids))))
  if (!length(ids)) return(.geo_empty())
  rows <- lapply(ids, function(gse) {
    if (is.na(gse) || !grepl("^GSE[0-9]+$", gse)) {
      return(.geo_row(gse, "invalid_accession", "Expected a GSE accession, e.g. GSE982."))
    }
    # GEO groups the last three digits: GSE982 -> GSEnnn; GSE1000 -> GSE1nnn.
    group <- sub("[0-9]{1,3}$", "nnn", gse)
    url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", group, "/", gse, "/matrix/")
    tryCatch({
      links <- .geo_retry(function() .geo_listing(url), retries, delay_seconds)
      files <- sort(unique(links[!is.na(links) & grepl(
        paste0("^", gse, "(-GPL[0-9]+)?_series_matrix[.]txt[.]gz$"), links
      )]))
      if (!length(files)) {
        return(.geo_row(gse, "no_matrix_files", "Listing contains no matching series matrices.", url = url))
      }
      do.call(rbind, lapply(files, function(file) {
        platform <- if (grepl("-GPL", file)) {
          sub(".*-(GPL[0-9]+)_.*", "\\1", file)
        } else NA_character_
        .geo_row(gse, "discovered", file = file, url = paste0(url, file), platform = platform)
      }))
    }, error = function(e) .geo_row(gse, "discovery_failed", conditionMessage(e), url = url),
    finally = Sys.sleep(delay_seconds))
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

# Stream the entire gzip to catch corruption without loading the matrix into RAM.
# A structurally valid matrix may still contain zero expression features.
.geo_validate_matrix <- function(path) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) {
    stop("Missing or empty file.")
  }
  magic <- readBin(path, what = "raw", n = 2L)
  if (!identical(magic, as.raw(c(31, 139)))) stop("File is not gzip-compressed.")
  con <- gzfile(path, "rt")
  on.exit(close(con))
  metadata <- begin <- end <- FALSE
  withCallingHandlers({
    repeat {
      lines <- readLines(con, n = 1000L, warn = TRUE)
      if (!length(lines)) break
      metadata <- metadata || any(grepl("^!(Series|Sample)_", lines))
      begin <- begin || any(grepl("^!series_matrix_table_begin", lines))
      end <- end || any(grepl("^!series_matrix_table_end", lines))
    }
  }, warning = function(w) stop(conditionMessage(w), call. = FALSE))
  if (!metadata || !begin || !end) stop("Missing series-matrix metadata or table boundaries.")
  invisible(TRUE)
}

.geo_accessions <- function(text, prefix = "GPL") {
  text <- as.character(unlist(text, use.names = FALSE))
  text <- text[!is.na(text)]
  unique(unlist(regmatches(text, gregexpr(paste0("\\b", prefix, "[0-9]+\\b"),
                                        text, perl = TRUE)), use.names = FALSE))
}

# Scan metadata up to the expression table, without a fixed line-count limit.
.geo_header_platforms <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  series <- samples <- gsm <- character()
  repeat {
    lines <- readLines(con, n = 100L, warn = FALSE)
    if (!length(lines)) break
    boundary <- which(grepl("^!series_matrix_table_begin", lines))
    if (length(boundary)) lines <- head(lines, boundary[1] - 1L)
    series <- c(series, grep("^!Series_platform_id[[:space:]=]", lines, value = TRUE))
    samples <- c(samples, grep("^!Sample_platform_id[[:space:]=]", lines, value = TRUE))
    gsm <- c(gsm, grep("^!Sample_geo_accession[[:space:]=]", lines, value = TRUE))
    if (length(boundary)) break
  }
  list(series = .geo_accessions(series), samples = .geo_accessions(samples),
       gsm = .geo_accessions(gsm, "GSM"))
}

# SOFT metadata avoids guessing which GPL belongs to a multi-platform matrix.
# Reference: https://seandavi.github.io/GEOquery/reference/getGEO.html
.geo_query_platform_metadata <- function(gse) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Platform fallback requires GEOquery. Install with BiocManager::install('GEOquery').")
  }
  object <- GEOquery::getGEO(gse, GSEMatrix = FALSE, getGPL = FALSE)
  samples <- GEOquery::GSMList(object)
  sample_platforms <- lapply(samples, function(sample) .geo_accessions(GEOquery::Meta(sample)$platform_id))
  list(series = .geo_accessions(c(GEOquery::Meta(object)$platform_id,
                                  unlist(sample_platforms, use.names = FALSE))),
       samples = sample_platforms)
}

# Also accepts an existing downloaded_matrices.csv read with read.csv().
# Download success and platform-lookup failure are reported separately.
fill_geo_platforms <- function(inventory, retries = 2L, delay_seconds = 0.4,
                               timeout_seconds = 300) {
  .geo_controls(retries, delay_seconds, timeout_seconds)
  stopifnot(is.data.frame(inventory),
            all(c("Platform", "SeriesAccession", "LocalFile", "Status") %in% names(inventory)))
  old <- options(timeout = max(getOption("timeout", 60), timeout_seconds))
  on.exit(options(old))
  result <- inventory
  result$Platform <- as.character(result$Platform)
  for (name in c("PlatformSource", "PlatformError")) {
    if (!name %in% names(result)) result[[name]] <- rep(NA_character_, nrow(result))
  }
  cache <- new.env(parent = emptyenv())
  for (i in seq_len(nrow(result))) {
    if (!is.na(result$Platform[i]) && nzchar(trimws(result$Platform[i]))) {
      if (is.na(result$PlatformSource[i]) || !nzchar(result$PlatformSource[i])) {
        result$PlatformSource[i] <- "existing_or_filename"
      }
      next
    }
    if (!result$Status[i] %in% c("downloaded", "already_exists")) next
    resolved <- tryCatch({
      path <- result$LocalFile[i]
      if (is.na(path) || !file.exists(path)) stop("Local matrix file is missing.")
      header <- .geo_header_platforms(path)
      if (length(header$samples)) {
        list(ids = header$samples, source = "matrix_sample_header")
      } else if (length(header$series) == 1L) {
        list(ids = header$series, source = "matrix_series_header")
      } else {
        gse <- result$SeriesAccession[i]
        if (is.na(gse) || !grepl("^GSE[0-9]+$", gse)) stop("Invalid GSE accession for platform lookup.")
        if (!exists(gse, envir = cache, inherits = FALSE)) {
          lookup <- tryCatch(list(value = .geo_retry(function() {
            on.exit(Sys.sleep(max(0.4, delay_seconds)))
            .geo_query_platform_metadata(gse)
          }, retries, delay_seconds)), error = function(e) list(error = conditionMessage(e)))
          assign(gse, lookup, envir = cache)
        }
        lookup <- get(gse, envir = cache, inherits = FALSE)
        if (!is.null(lookup$error)) stop(lookup$error)
        meta <- lookup$value
        matched <- meta$samples[header$gsm]
        if (length(header$gsm) && all(lengths(matched) > 0L)) {
          list(ids = unique(unlist(matched, use.names = FALSE)), source = "GEOquery_samples")
        } else if (length(meta$series) == 1L) {
          list(ids = meta$series, source = "GEOquery_series")
        } else {
          stop("No unambiguous platform mapping: GEO returned zero or multiple platforms without a complete matrix sample match.")
        }
      }
    }, error = function(e) list(error = conditionMessage(e)))
    result$PlatformError[i] <- if (is.null(resolved$error)) NA_character_ else resolved$error
    result$PlatformSource[i] <- if (is.null(resolved$error)) resolved$source else "unresolved"
    if (is.null(resolved$error)) result$Platform[i] <- paste(resolved$ids, collapse = "; ")
  }
  result
}

# Retains discovery failures in the returned inventory. Each file is independent.
# Valid existing files are reused; invalid existing files are downloaded again.
# overwrite=TRUE refreshes even valid files, retaining the old file if fetching fails.
download_geo_matrices <- function(inventory, dest_dir, overwrite = FALSE,
                                  retries = 2L, delay_seconds = 0.4,
                                  timeout_seconds = 300, resolve_platforms = TRUE) {
  .geo_controls(retries, delay_seconds, timeout_seconds)
  if (!is.data.frame(inventory) || !all(names(.geo_empty()) %in% names(inventory))) {
    stop("inventory must be returned by discover_geo_matrices().")
  }
  stopifnot(length(overwrite) == 1L, !is.na(overwrite), is.logical(overwrite))
  stopifnot(length(resolve_platforms) == 1L, !is.na(resolve_platforms), is.logical(resolve_platforms))
  if (!dir.exists(dest_dir) && !dir.create(dest_dir, recursive = TRUE)) {
    stop("Cannot create destination directory: ", dest_dir)
  }
  dest_dir <- normalizePath(dest_dir, mustWork = TRUE)
  old <- options(timeout = max(getOption("timeout", 60), timeout_seconds))
  on.exit(options(old))
  result <- inventory
  eligible <- c("discovered", "downloaded", "already_exists", "download_failed")
  for (i in seq_len(nrow(result))) {
    if (is.na(result$Status[i]) || !result$Status[i] %in% eligible) next
    row <- result[i, ]
    result$LocalFile[i] <- NA_character_
    outcome <- tryCatch({
      if (is.na(row$MatrixFile) || !grepl(
        "^GSE[0-9]+(-GPL[0-9]+)?_series_matrix[.]txt[.]gz$", row$MatrixFile
      ) || is.na(row$URL)) stop("Invalid matrix filename or missing URL.")
      target <- file.path(dest_dir, row$MatrixFile)
      existing_valid <- file.exists(target) && isTRUE(tryCatch({
        .geo_validate_matrix(target)
        TRUE
      }, error = function(e) FALSE))
      if (existing_valid && !overwrite) {
        list(status = "already_exists", path = target)
      } else {
        .geo_retry(function() {
          # Same directory permits an atomic rename on supported filesystems.
          partial <- tempfile(pattern = paste0(".", row$MatrixFile, "."), tmpdir = dest_dir)
          on.exit(unlink(partial))
          .geo_fetch(row$URL, partial)
          .geo_validate_matrix(partial)
          if (!file.rename(partial, target)) stop("Cannot move completed download to: ", target)
        }, retries, delay_seconds)
        list(status = "downloaded", path = target)
      }
    }, error = function(e) list(status = "download_failed", error = conditionMessage(e)),
    finally = Sys.sleep(delay_seconds))
    result$Status[i] <- outcome$status
    result$ErrorMessage[i] <- if (is.null(outcome$error)) NA_character_ else outcome$error
    if (!is.null(outcome$path)) result$LocalFile[i] <- outcome$path
  }
  if (resolve_platforms) {
    result <- fill_geo_platforms(result, retries = retries, delay_seconds = delay_seconds,
                                timeout_seconds = timeout_seconds)
  }
  result
}

# Diagnosis abbreviation to full name, shared with 01_get_AML_GEO_candidates.R.
GEO_DIAGNOSIS_MAP <- c(
  AML = "acute myeloid leukemia",
  ALL = "acute lymphocytic leukemia"
)

# Build the DMSO/vehicle query strings for a supported diagnosis abbreviation.
geo_diagnosis_query_strings <- function(abbreviation, diagnosis_map = GEO_DIAGNOSIS_MAP) {
  if (!abbreviation %in% names(diagnosis_map)) {
    stop(
      "Unsupported diagnosis '", abbreviation, "'. Supported values: ",
      paste(names(diagnosis_map), collapse = ", ")
    )
  }
  full_name <- diagnosis_map[[abbreviation]]
  c(
    paste(abbreviation, "DMSO"),
    paste0("\"", full_name, "\" DMSO"),
    paste(abbreviation, "vehicle"),
    paste0("\"", full_name, "\" vehicle")
  )
}

# Build a query from literal words/phrases; use query_strings directly for complex
# Entrez expressions, e.g. '(AML OR "acute myeloid leukemia") AND (DMSO OR vehicle)'.
geo_keyword_query <- function(keywords, match = c("all", "any"), organism = NULL) {
  match <- match.arg(match)
  keywords <- unique(trimws(as.character(keywords)))
  if (!length(keywords) || anyNA(keywords) || any(!nzchar(keywords)) ||
      any(grepl('"', keywords, fixed = TRUE))) stop("Provide nonempty keywords without embedded quotes.")
  query <- paste0("(", paste(paste0('"', keywords, '"[All Fields]'),
                             collapse = if (match == "all") " AND " else " OR "), ")")
  if (!is.null(organism)) {
    stopifnot(length(organism) == 1L, !is.na(organism), nzchar(organism),
              !grepl('"', organism, fixed = TRUE))
    query <- paste0(query, ' AND "', organism, '"[Organism]')
  }
  query
}

.geo_eutils <- function(endpoint, params, retries, delay_seconds) {
  url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", endpoint, ".fcgi?",
                paste(paste0(names(params), "=", vapply(params, function(x) {
                  utils::URLencode(as.character(x), reserved = TRUE)
                }, character(1))), collapse = "&"))
  .geo_retry(function() {
    file <- tempfile(fileext = ".xml")
    on.exit(unlink(file))
    on.exit(Sys.sleep(delay_seconds), add = TRUE)
    .geo_fetch(url, file)
    doc <- xml2::read_xml(file)
    errors <- xml2::xml_text(xml2::xml_find_all(doc, "//ERROR | //ErrorList/*"))
    if (length(errors)) stop(paste(errors, collapse = "; "))
    doc
  }, retries, delay_seconds)
}

.geo_xml_text <- function(doc, xpath) {
  xml2::xml_text(xml2::xml_find_first(doc, xpath), trim = TRUE)
}

.geo_study_row <- function(uid, accession = NA_character_, title = NA_character_,
                           summary = NA_character_, organism = NA_character_,
                           type = NA_character_, status = "metadata_failed",
                           error = NA_character_) {
  data.frame(UID = uid, SeriesAccession = accession, Title = title, Summary = summary,
             Organism = organism, Type = type, MetadataStatus = status,
             MetadataError = error, stringsAsFactors = FALSE)
}

# Each query is executed separately and overlapping studies are deduplicated.
# keywords are literal, case-insensitive substring annotations of title/summary;
# they do not filter the Entrez search results or prove sample-level treatment.
# A finite max_studies_per_query is explicitly reported as 'truncated'.
discover_geo_by_keywords <- function(query_strings, keywords = character(),
                                     max_studies_per_query = Inf, page_size = 200L,
                                     retries = 2L, delay_seconds = 0.4,
                                     timeout_seconds = 120) {
  if (!exists("discover_geo_matrices", mode = "function")) {
    stop("Matrix discovery functions are unavailable.")
  }
  if (!requireNamespace("xml2", quietly = TRUE)) stop("Install xml2 first.")
  .geo_controls(retries, delay_seconds, timeout_seconds)
  # This implementation uses unauthenticated E-utilities requests.
  delay_seconds <- max(0.4, delay_seconds)
  stopifnot(length(max_studies_per_query) == 1L, !is.na(max_studies_per_query),
            max_studies_per_query > 0,
            max_studies_per_query == floor(max_studies_per_query),
            length(page_size) == 1L, is.finite(page_size), page_size >= 1,
            page_size <= 500L, page_size == floor(page_size))
  queries <- unique(trimws(as.character(query_strings)))
  if (!length(queries) || anyNA(queries) || any(!nzchar(queries))) stop("Provide nonempty search queries.")
  keywords <- unique(trimws(as.character(keywords)))
  keywords <- keywords[!is.na(keywords) & nzchar(keywords)]
  old <- options(timeout = max(getOption("timeout", 60), timeout_seconds))
  on.exit(options(old))
  request <- function(endpoint, params) .geo_eutils(endpoint, params, retries, delay_seconds)
  links <- data.frame(QueryString = character(), UID = character())
  logs <- vector("list", length(queries))
  for (i in seq_along(queries)) {
    query <- paste0("(", queries[i], ") AND gse[ETYP]")
    ids <- character()
    total <- NA_real_
    translated <- NA_character_
    warnings <- character()
    error <- NA_character_
    state <- tryCatch({
      offset <- 0L
      repeat {
        n <- min(page_size, max_studies_per_query - offset)
        doc <- request("esearch", list(db = "gds", term = query, retmode = "xml",
                                        retstart = offset, retmax = n))
        count <- suppressWarnings(as.numeric(.geo_xml_text(doc, "/eSearchResult/Count")))
        if (is.na(count)) stop("Missing ESearch result count.")
        if (is.na(total)) total <- count
        translated <- .geo_xml_text(doc, "/eSearchResult/QueryTranslation")
        warnings <- unique(c(warnings, xml2::xml_text(xml2::xml_find_all(doc, "//WarningList/*"))))
        page <- xml2::xml_text(xml2::xml_find_all(doc, "/eSearchResult/IdList/Id"))
        ids <- unique(c(ids, page))
        offset <- offset + length(page)
        if (offset >= min(total, max_studies_per_query)) break
        if (!length(page)) stop("ESearch pagination stopped before all requested records were returned.")
      }
      if (total == 0) "no_hits" else if (total > max_studies_per_query) "truncated" else "complete"
    }, error = function(e) {
      error <<- conditionMessage(e)
      if (length(ids)) "partial_failure" else "search_failed"
    })
    if (length(ids)) links <- rbind(links, data.frame(QueryString = queries[i], UID = ids))
    logs[[i]] <- data.frame(QueryString = queries[i], EntrezQuery = query,
                            QueryTranslation = translated, TotalHits = total,
                            RetrievedUIDs = length(ids), Status = state,
                            Warnings = paste(warnings, collapse = "; "), ErrorMessage = error)
  }
  uids <- unique(links$UID)
  studies <- .geo_study_row(NA_character_)[FALSE, ]
  batches <- split(uids, ceiling(seq_along(uids) / 100L))
  for (batch in batches) {
    batch_error <- NA_character_
    doc <- tryCatch(request("esummary", list(db = "gds", id = paste(batch, collapse = ","),
                                              retmode = "xml")), error = function(e) {
      batch_error <<- conditionMessage(e)
      NULL
    })
    rows <- lapply(batch, function(uid) {
      if (is.null(doc)) return(.geo_study_row(uid, error = batch_error))
      if (!grepl("^[0-9]+$", uid)) return(.geo_study_row(uid, error = "Invalid Entrez UID."))
      record <- xml2::xml_find_first(doc, paste0("//DocSum[Id='", uid, "']"))
      if (inherits(record, "xml_missing")) return(.geo_study_row(uid, error = "Summary missing from response."))
      item <- function(name) .geo_xml_text(record, paste0(".//Item[@Name='", name, "']"))
      acc <- item("Accession")
      if (is.na(acc) || !grepl("^GSE[0-9]+$", acc)) {
        return(.geo_study_row(uid, error = "Summary has no valid GSE accession."))
      }
      .geo_study_row(uid, acc, item("title"), item("summary"), item("taxon"),
                     item("gdsType"), status = "OK")
    })
    studies <- rbind(studies, do.call(rbind, rows))
  }
  studies$QueryStrings <- vapply(studies$UID, function(uid) {
    paste(unique(links$QueryString[links$UID == uid]), collapse = "; ")
  }, character(1))
  hits <- function(text) {
    if (is.na(text)) return(NA_character_)
    found <- keywords[vapply(keywords, function(k) grepl(tolower(k), tolower(text),
                                                       fixed = TRUE), logical(1))]
    if (length(found)) paste(found, collapse = "; ") else NA_character_
  }
  studies$TitleKeywordHits <- vapply(studies$Title, hits, character(1))
  studies$SummaryKeywordHits <- vapply(studies$Summary, hits, character(1))
  valid <- studies$MetadataStatus == "OK"
  matrices <- discover_geo_matrices(studies$SeriesAccession[valid], retries = retries,
                                    delay_seconds = delay_seconds, timeout_seconds = timeout_seconds)
  annotation <- match(matrices$SeriesAccession, studies$SeriesAccession)
  for (field in c("Title", "Summary", "QueryStrings", "TitleKeywordHits", "SummaryKeywordHits")) {
    matrices[[field]] <- studies[[field]][annotation]
  }
  list(studies = studies, matrices = matrices, query_log = do.call(rbind, logs),
       query_study_links = links)
}

.geo_bind <- function(a, b) {
  columns <- union(names(a), names(b))
  for (name in setdiff(columns, names(a))) a[[name]] <- rep(NA, nrow(a))
  for (name in setdiff(columns, names(b))) b[[name]] <- rep(NA, nrow(b))
  rbind(a[columns], b[columns])
}

.geo_atomic_csv <- function(data, path) {
  temporary <- tempfile("report-", tmpdir = dirname(path), fileext = ".csv")
  on.exit(unlink(temporary))
  write.csv(data, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) stop("Cannot replace report: ", path)
}

# Extract actual matrix sample counts and study metadata, without reading values.
.geo_report_metadata <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  header <- character()
  repeat {
    lines <- readLines(con, n = 100L, warn = FALSE)
    if (!length(lines)) break
    boundary <- which(grepl("^!series_matrix_table_begin", lines))
    if (length(boundary)) lines <- head(lines, boundary[1] - 1L)
    header <- c(header, lines[grepl("^!(Series_title|Sample_organism_ch1|Sample_geo_accession)([[:space:]]|=)", lines)])
    if (length(boundary)) break
  }
  values <- function(key) {
    lines <- grep(paste0("^!", key, "([[:space:]]|=)"), header, value = TRUE)
    text <- sub(paste0("^!", key, "[[:space:]=]+"), "", lines)
    unique(trimws(gsub('"', "", unlist(strsplit(text, "\t", fixed = TRUE)), fixed = TRUE)))
  }
  title <- values("Series_title")
  organism <- values("Sample_organism_ch1")
  gsm <- .geo_accessions(values("Sample_geo_accession"), "GSM")
  list(SeriesTitle = if (length(title)) paste(title, collapse = "; ") else NA_character_,
       Organism = if (length(organism)) paste(organism, collapse = "; ") else NA_character_,
       SampleCount = if (length(gsm)) length(gsm) else NA_integer_)
}

.geo_report_key <- function(data) {
  paste(data$SeriesAccession, ifelse(is.na(data$MatrixFile) | data$MatrixFile == "", "__study__", data$MatrixFile), sep = "|")
}

# Reports merge by study/file, retaining records for unprocessed studies on reruns.
# Every validated file is checkpointed, not just each completed study.
run_geo_pipeline <- function(input = "GSE_Metadata_Inventory.xlsx",
                             out_dir = "downloads/geo_metadata_inventory",
                             query_strings = NULL, diagnosis = NULL, keywords = character(),
                             max_studies_per_query = Inf, retries = 1L) {
  if (!is.null(query_strings) && !is.null(diagnosis)) {
    stop("Provide either query_strings or diagnosis, not both.")
  }
  if (is.null(query_strings) && !is.null(diagnosis)) {
    query_strings <- geo_diagnosis_query_strings(diagnosis)
    if (!length(keywords)) {
      keywords <- c(diagnosis, GEO_DIAGNOSIS_MAP[[diagnosis]], "DMSO", "vehicle")
    }
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = TRUE)
  matrix_dir <- file.path(out_dir, "matrices")
  dir.create(matrix_dir, showWarnings = FALSE)
  started <- Sys.time()
  prior_path <- file.path(out_dir, "all_results.csv")
  reports <- if (file.exists(prior_path)) {
    read.csv(prior_path, stringsAsFactors = FALSE, colClasses = "character")
  } else .geo_empty()
  stopifnot(all(names(.geo_empty()) %in% names(reports)))
  # Also preserve any additional rows from a previously maintained success CSV.
  success_path <- file.path(out_dir, "downloaded_matrices.csv")
  if (file.exists(success_path)) {
    prior_success <- read.csv(success_path, stringsAsFactors = FALSE, colClasses = "character")
    stopifnot(all(names(.geo_empty()) %in% names(prior_success)))
    reports <- .geo_bind(reports[!.geo_report_key(reports) %in% .geo_report_key(prior_success), ], prior_success)
  }
  reports <- reports[!duplicated(.geo_report_key(reports), fromLast = TRUE), ]
  discovered <- NULL
  if (is.null(query_strings)) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl for Excel input.")
    sheets <- readxl::excel_sheets(input)
    # Accept the inventory workbook produced by heme_cancer_geo.py. It stores
    # study/sample rows across cancer tabs and calls the accession column GSE.
    sheet <- if ("Metadata" %in% sheets) "Metadata" else NULL
    if (is.null(sheet)) {
      usable <- setdiff(sheets, c("README", "Study_Index", "Export_Notes"))
      if (!length(usable)) stop("Workbook has no metadata sheets.")
      parts <- lapply(usable, function(s) as.data.frame(readxl::read_excel(input, sheet = s)))
      metadata <- do.call(rbind, lapply(parts, function(x) {
        if ("GSE" %in% names(x) && !("SeriesAccession" %in% names(x))) x$SeriesAccession <- x$GSE
        if ("StudyTitle" %in% names(x) && !("SeriesTitle" %in% names(x))) x$SeriesTitle <- x$StudyTitle
        if ("PubMedIDs" %in% names(x) && !("PubMedID" %in% names(x))) x$PubMedID <- x$PubMedIDs
        if ("SeriesSampleCount" %in% names(x) && !("SampleCount" %in% names(x))) x$SampleCount <- x$SeriesSampleCount
        x
      }))
    } else metadata <- as.data.frame(readxl::read_excel(input, sheet = sheet))
    source_name <- basename(input)
  } else {
    search <- discover_geo_by_keywords(query_strings, keywords, max_studies_per_query,
                                       retries = retries)
    for (name in names(search)) .geo_atomic_csv(search[[name]], file.path(out_dir, paste0("search_", name, ".csv")))
    metadata <- search$studies[search$studies$MetadataStatus == "OK", ]
    metadata$SeriesTitle <- metadata$Title
    discovered <- search$matrices
    source_name <- "GEO keyword search"
  }
  stopifnot("SeriesAccession" %in% names(metadata))
  metadata$SeriesAccession <- toupper(trimws(as.character(metadata$SeriesAccession)))
  ids <- unique(metadata$SeriesAccession)
  processed <- 0L
  save_reports <- function() {
    ok <- reports$Status %in% c("downloaded", "already_exists")
    missing_file <- ok & (is.na(reports$LocalFile) | !file.exists(reports$LocalFile))
    reports$Status[missing_file] <<- "local_file_missing"
    reports$ErrorMessage[missing_file] <<- "Previously recorded matrix is missing; rerun to download again."
    ok <- reports$Status %in% c("downloaded", "already_exists")
    blank <- is.na(reports$Platform) | !nzchar(trimws(reports$Platform))
    reports$FileSizeBytes <<- rep(NA_real_, nrow(reports))
    reports$FileSizeBytes[ok] <<- file.info(reports$LocalFile[ok])$size
    .geo_atomic_csv(reports, prior_path)
    .geo_atomic_csv(reports[ok, ], success_path)
    .geo_atomic_csv(reports[!ok, ], file.path(out_dir, "failures.csv"))
    .geo_atomic_csv(reports[ok & blank, ], file.path(out_dir, "platform_failures.csv"))
    summary <- data.frame(
      Metric = c("Unique requested accessions", "Accessions processed this run", "Report rows (including prior runs)",
                 "Validated local matrix files", "Files with platform IDs", "Files missing platform IDs",
                 "Failure rows", "Downloaded bytes", "Elapsed seconds this run"),
      Value = c(length(ids), processed, nrow(reports), sum(ok), sum(ok & !blank), sum(ok & blank),
                sum(!ok), sum(reports$FileSizeBytes[ok]), as.numeric(difftime(Sys.time(), started, units = "secs"))))
    .geo_atomic_csv(summary, file.path(out_dir, "summary.csv"))
  }
  checkpoint <- function(rows) {
    # Clear obsolete study-level discovery errors once a matrix is available.
    new_files <- rows$SeriesAccession[!is.na(rows$MatrixFile) & nzchar(rows$MatrixFile)]
    stale_error <- reports$SeriesAccession %in% new_files &
      (is.na(reports$MatrixFile) | !nzchar(reports$MatrixFile))
    reports <<- .geo_bind(reports[!stale_error & !.geo_report_key(reports) %in% .geo_report_key(rows), ], rows)
    save_reports()
  }
  for (i in seq_along(ids)) {
    message("[", i, "/", length(ids), "] ", ids[i])
    found <- if (is.null(discovered)) discover_geo_matrices(ids[i], retries = retries) else {
      discovered[discovered$SeriesAccession == ids[i], ]
    }
    for (j in seq_len(nrow(found))) {
      row <- download_geo_matrices(found[j, ], matrix_dir, retries = retries)
      index <- match(row$SeriesAccession, metadata$SeriesAccession)
      for (field in intersect(c("SeriesTitle", "Organism", "PubMedID", "SampleCount", "QueryStrings",
                                "TitleKeywordHits", "SummaryKeywordHits"), names(metadata))) {
        row[[field]] <- metadata[[field]][index]
      }
      row$MetadataError <- NA_character_
      if (row$Status %in% c("downloaded", "already_exists")) {
        details <- tryCatch(.geo_report_metadata(row$LocalFile), error = function(e) {
          row$MetadataError <<- conditionMessage(e)
          list()
        })
        for (field in names(details)) if (!is.na(details[[field]])) row[[field]] <- details[[field]]
      }
      row$SourceInput <- source_name
      row$ProcessedAtUTC <- format(Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ")
      checkpoint(row)
    }
    processed <- i
    save_reports()
  }
  save_reports()
  message("Complete. Reports: ", out_dir)
  invisible(reports)
}

.geo_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) && args[1] == "--help") {
    cat("Workbook: Rscript 05_Download_Metadata_Inventory.R [input.xlsx] [output_directory]\n",
        "Search:   Rscript 05_Download_Metadata_Inventory.R --query '<condition> AND DMSO' --out downloads/geo_search\n",
        paste0("Diagnosis: Rscript 05_Download_Metadata_Inventory.R --diagnosis <",
               paste(names(GEO_DIAGNOSIS_MAP), collapse = "|"), "> --out downloads/geo_search\n"),
        "Optional search arguments: --keywords '<term1>,<term2>' --max-studies 10\n")
    return(invisible(NULL))
  }
  if (length(args) && args[1] %in% c("--query", "--diagnosis")) {
    if (length(args) %% 2L != 0L || any(!args[seq(1L, length(args), 2L)] %in%
                                        c("--query", "--diagnosis", "--out", "--keywords", "--max-studies"))) stop("Invalid arguments; use --help.")
    options <- setNames(as.list(args[seq(2L, length(args), 2L)]), args[seq(1L, length(args), 2L)])
    get_option <- function(key, default) if (is.null(options[[key]])) default else options[[key]]
    keyword_option <- get_option("--keywords", "")
    run_geo_pipeline(query_strings = options[["--query"]],
                     diagnosis = options[["--diagnosis"]],
                     out_dir = get_option("--out", "downloads/geo_search"),
                     keywords = if (nzchar(keyword_option)) strsplit(keyword_option, ",", fixed = TRUE)[[1]] else character(),
                     max_studies_per_query = as.numeric(get_option("--max-studies", "Inf")))
  } else {
    if (length(args) > 2L) stop("Invalid arguments; use --help.")
    run_geo_pipeline(input = if (length(args)) args[1] else "GSE_Metadata_Inventory.xlsx",
                     out_dir = if (length(args) > 1L) args[2] else "downloads/geo_metadata_inventory")
  }
}

if (sys.nframe() == 0L) .geo_main()
