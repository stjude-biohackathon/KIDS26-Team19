## Data helpers for the Heme GEO download planner Shiny app.
## Reads Heme_Cancer_GEO_Inventory.xlsx, filters studies, and reports bronze downloads.

HEME_SAMPLE_SHEETS <- c(
  "Pediatric_AML",
  "Adult_AML",
  "ALL",
  "CLL",
  "Review_Unassigned"
)

HEME_REQUIRED_COLUMNS <- c(
  "CancerType",
  "GSE",
  "GSM",
  "StudyTitle",
  "Organism",
  "SpecimenType",
  "SampleType",
  "TreatmentStatus",
  "TreatmentAgents",
  "CellLineNames",
  "SeriesMatrixURL",
  "SeriesSampleCount"
)

HEME_BLANK_TOKENS <- c(
  "",
  "not reported",
  "none detected",
  "unknown",
  "na",
  "n/a"
)

heme_blank_like <- function(x) {
  text <- trimws(as.character(x))
  is.na(x) | !nzchar(text) | tolower(text) %in% HEME_BLANK_TOKENS
}

heme_split_tokens <- function(x) {
  text <- as.character(x)
  text[is.na(text)] <- ""
  parts <- unlist(strsplit(text, "\\s*\\|\\s*"), use.names = FALSE)
  parts <- trimws(parts)
  parts[!heme_blank_like(parts)]
}

heme_normalize_cell_line <- function(x) {
  toupper(gsub("[^A-Za-z0-9]+", "", as.character(x)))
}

heme_unique_sorted <- function(x) {
  values <- unique(as.character(x))
  values <- values[!heme_blank_like(values)]
  sort(values)
}

default_heme_inventory_path <- function(repo_root) {
  candidates <- c(
    file.path(repo_root, "R_Scripts", "Heme_Cancer_GEO_Inventory.xlsx"),
    file.path(dirname(repo_root), "KIDS26-Team19-main", "R_Scripts", "Heme_Cancer_GEO_Inventory.xlsx")
  )
  found <- candidates[file.exists(candidates)]
  if (length(found)) found[[1]] else candidates[[1]]
}

read_heme_inventory <- function(path) {
  if (!file.exists(path)) {
    stop("Heme inventory not found: ", path)
  }
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Install the readxl package first.")
  }
  sheets <- readxl::excel_sheets(path)
  sample_sheets <- intersect(HEME_SAMPLE_SHEETS, sheets)
  if (!length(sample_sheets)) {
    stop(
      "No sample sheets found. Expected one of: ",
      paste(HEME_SAMPLE_SHEETS, collapse = ", ")
    )
  }
  rows <- lapply(sample_sheets, function(sheet) {
    data <- as.data.frame(readxl::read_excel(path, sheet = sheet))
    missing <- setdiff(HEME_REQUIRED_COLUMNS, names(data))
    if (length(missing)) {
      stop("Sheet ", sheet, " is missing columns: ", paste(missing, collapse = ", "))
    }
    if (!nrow(data)) {
      return(NULL)
    }
    numeric_cols <- intersect(c("SeriesSampleCount", "AgeYears"), names(data))
    for (col in setdiff(names(data), numeric_cols)) {
      data[[col]] <- as.character(data[[col]])
    }
    for (col in numeric_cols) {
      data[[col]] <- suppressWarnings(as.numeric(data[[col]]))
    }
    data$InventorySheet <- sheet
    data$GSE <- toupper(trimws(as.character(data$GSE)))
    data
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    stop("No sample rows found in ", path)
  }
  inventory <- dplyr::bind_rows(rows)
  inventory <- as.data.frame(inventory, stringsAsFactors = FALSE)
  rownames(inventory) <- NULL
  inventory
}

heme_filter_choices <- function(inventory) {
  cell_lines <- heme_unique_sorted(heme_split_tokens(inventory$CellLineNames))
  list(
    diseases = heme_unique_sorted(inventory$CancerType),
    cell_lines = cell_lines,
    organisms = heme_unique_sorted(inventory$Organism),
    specimen_types = heme_unique_sorted(inventory$SpecimenType),
    treatment_statuses = heme_unique_sorted(inventory$TreatmentStatus)
  )
}

.heme_selected <- function(selected) {
  if (is.null(selected)) {
    return(character())
  }
  values <- unique(trimws(as.character(selected)))
  values[!heme_blank_like(values)]
}

.heme_matches_any <- function(values, selected) {
  selected <- .heme_selected(selected)
  if (!length(selected)) {
    return(rep(TRUE, length(values)))
  }
  as.character(values) %in% selected
}

.heme_row_cell_lines_match <- function(cell_line_names, selected) {
  selected <- .heme_selected(selected)
  if (!length(selected)) {
    return(rep(TRUE, length(cell_line_names)))
  }
  wanted <- unique(heme_normalize_cell_line(selected))
  wanted <- wanted[nzchar(wanted)]
  vapply(as.character(cell_line_names), function(value) {
    tokens <- heme_normalize_cell_line(heme_split_tokens(value))
    any(tokens %in% wanted)
  }, logical(1))
}

.heme_vehicle_series <- function(inventory) {
  sample_type <- as.character(inventory$SampleType)
  status <- as.character(inventory$TreatmentStatus)
  agents <- as.character(inventory$TreatmentAgents)
  has_vehicle <- grepl("DMSO", sample_type, ignore.case = TRUE) |
    grepl("vehicle", status, ignore.case = TRUE) |
    grepl("DMSO|vehicle", agents, ignore.case = TRUE)
  unique(inventory$GSE[has_vehicle %in% TRUE])
}

filter_heme_inventory <- function(inventory,
                                  diseases = NULL,
                                  cell_lines = NULL,
                                  organisms = NULL,
                                  specimen_types = NULL,
                                  treatment_statuses = NULL,
                                  require_vehicle_control = FALSE) {
  keep <- .heme_matches_any(inventory$CancerType, diseases) &
    .heme_row_cell_lines_match(inventory$CellLineNames, cell_lines) &
    .heme_matches_any(inventory$Organism, organisms) &
    .heme_matches_any(inventory$SpecimenType, specimen_types) &
    .heme_matches_any(inventory$TreatmentStatus, treatment_statuses)
  filtered <- inventory[keep, , drop = FALSE]
  if (isTRUE(require_vehicle_control)) {
    vehicle_gse <- .heme_vehicle_series(inventory)
    filtered <- filtered[filtered$GSE %in% vehicle_gse, , drop = FALSE]
  }
  rownames(filtered) <- NULL
  filtered
}

heme_download_plan <- function(samples) {
  if (!nrow(samples)) {
    return(data.frame(
      GSE = character(),
      Diseases = character(),
      StudyTitle = character(),
      Organism = character(),
      CellLines = character(),
      MatchingSamples = integer(),
      SeriesSampleCount = integer(),
      MatrixFileCount = integer(),
      SeriesMatrixURL = character(),
      stringsAsFactors = FALSE
    ))
  }
  split_urls <- function(url) {
    unique(heme_split_tokens(url))
  }
  gse <- unique(samples$GSE)
  rows <- lapply(gse, function(id) {
    part <- samples[samples$GSE == id, , drop = FALSE]
    urls <- unique(unlist(lapply(part$SeriesMatrixURL, split_urls), use.names = FALSE))
    data.frame(
      GSE = id,
      Diseases = paste(sort(unique(part$CancerType)), collapse = "; "),
      StudyTitle = paste(unique(part$StudyTitle[!heme_blank_like(part$StudyTitle)]), collapse = "; "),
      Organism = paste(sort(unique(part$Organism[!heme_blank_like(part$Organism)])), collapse = "; "),
      CellLines = paste(heme_unique_sorted(heme_split_tokens(part$CellLineNames)), collapse = "; "),
      MatchingSamples = nrow(part),
      SeriesSampleCount = suppressWarnings(max(as.numeric(part$SeriesSampleCount), na.rm = TRUE)),
      MatrixFileCount = length(urls),
      SeriesMatrixURL = paste(urls, collapse = " | "),
      stringsAsFactors = FALSE
    )
  })
  plan <- do.call(rbind, rows)
  plan$SeriesSampleCount[!is.finite(plan$SeriesSampleCount)] <- NA_real_
  plan[order(plan$GSE), , drop = FALSE]
}

heme_plan_files <- function(plan) {
  if (!nrow(plan)) {
    return(data.frame(
      GSE = character(),
      URL = character(),
      MatrixFile = character(),
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(seq_len(nrow(plan)), function(i) {
    urls <- heme_split_tokens(plan$SeriesMatrixURL[i])
    data.frame(
      GSE = rep(plan$GSE[i], length(urls)),
      URL = urls,
      MatrixFile = basename(urls),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

empty_bronze_manifest <- function() {
  data.frame(
    SeriesAccession = character(),
    Platform = character(),
    MatrixFile = character(),
    URL = character(),
    Bytes = numeric(),
    MD5 = character(),
    DownloadedAt = character(),
    Status = character(),
    stringsAsFactors = FALSE
  )
}

read_download_report <- function(bronze_dir = "bronze") {
  path <- file.path(bronze_dir, "manifest.csv")
  if (!file.exists(path)) {
    report <- empty_bronze_manifest()
    report$LocalFile <- character()
    report$LocalFileExists <- logical()
    return(report)
  }
  report <- read.csv(path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  required <- names(empty_bronze_manifest())
  missing <- setdiff(required, names(report))
  if (length(missing)) {
    stop("bronze/manifest.csv is missing columns: ", paste(missing, collapse = ", "))
  }
  report$LocalFile <- ifelse(
    heme_blank_like(report$MatrixFile),
    NA_character_,
    file.path(bronze_dir, "geo", report$SeriesAccession, report$MatrixFile)
  )
  report$LocalFileExists <- !is.na(report$LocalFile) & file.exists(report$LocalFile)
  report
}

summarize_download_report <- function(report) {
  ok <- !is.na(report$Status) & report$Status == "OK"
  data.frame(
    Metric = c(
      "Manifest rows",
      "Unique studies",
      "Successful files",
      "Failed files",
      "Local files still on disk",
      "Downloaded bytes"
    ),
    Value = c(
      nrow(report),
      length(unique(report$SeriesAccession[nzchar(as.character(report$SeriesAccession))])),
      sum(ok),
      sum(!ok),
      sum(report$LocalFileExists %in% TRUE),
      sum(as.numeric(report$Bytes[ok]), na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}

merge_plan_with_downloads <- function(plan, manifest) {
  if (!nrow(plan)) {
    plan$AlreadyDownloaded <- logical()
    plan$DownloadedFiles <- integer()
    plan$DownloadStatus <- character()
    return(plan)
  }
  if (!nrow(manifest)) {
    plan$AlreadyDownloaded <- FALSE
    plan$DownloadedFiles <- 0L
    plan$DownloadStatus <- "not downloaded"
    return(plan)
  }
  ok <- !is.na(manifest$Status) &
    manifest$Status %in% c("OK", "downloaded", "already_exists")
  counts <- tapply(ok, manifest$SeriesAccession, sum)
  status_by_gse <- tapply(as.character(manifest$Status), manifest$SeriesAccession, function(x) {
    paste(sort(unique(x)), collapse = "; ")
  })
  plan$DownloadedFiles <- as.integer(unname(counts[plan$GSE]))
  plan$DownloadedFiles[is.na(plan$DownloadedFiles)] <- 0L
  plan$AlreadyDownloaded <- plan$DownloadedFiles > 0L
  plan$DownloadStatus <- unname(as.character(status_by_gse[plan$GSE]))
  plan$DownloadStatus[is.na(plan$DownloadStatus)] <- "not downloaded"
  plan
}

heme_run_downloads <- function(gse_ids,
                               bronze_dir = "bronze",
                               delay_seconds = 0.4,
                               download_fun = NULL,
                               progress_fun = NULL) {
  ids <- unique(toupper(trimws(as.character(gse_ids))))
  ids <- ids[!heme_blank_like(ids)]
  if (!length(ids)) {
    stop("No GSE accessions selected.")
  }
  if (is.null(download_fun)) {
    if (!exists("download_study", mode = "function")) {
      stop("download_study() is not available. Source R/bronze.R first.")
    }
    download_fun <- download_study
  }
  dir.create(bronze_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    if (is.function(progress_fun)) {
      progress_fun(i, length(ids), ids[i])
    }
    rows[[i]] <- download_fun(ids[i], bronze.dir = bronze_dir)
    if (i < length(ids) && delay_seconds > 0) {
      Sys.sleep(delay_seconds)
    }
  }
  new_rows <- do.call(rbind, rows)
  manifest_path <- file.path(bronze_dir, "manifest.csv")
  if (file.exists(manifest_path)) {
    prior <- read.csv(manifest_path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
    prior <- prior[!prior$SeriesAccession %in% ids, , drop = FALSE]
    manifest <- rbind(prior, new_rows)
  } else {
    manifest <- new_rows
  }
  write.csv(manifest, manifest_path, row.names = FALSE, na = "")
  manifest
}

heme_plan_summary_text <- function(plan) {
  if (!nrow(plan)) {
    return("No studies match the current filters.")
  }
  paste0(
    nrow(plan), " studies, ",
    sum(plan$MatchingSamples), " matching samples, ",
    sum(plan$MatrixFileCount), " series-matrix files."
  )
}

PIPELINE_OK_STATUS <- c("OK", "downloaded", "already_exists")

read_csv_if_exists <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  read.csv(path, stringsAsFactors = FALSE, na.strings = c("", "NA"), check.names = FALSE)
}

read_pipeline_manifest <- function(out_dir) {
  success <- read_csv_if_exists(file.path(out_dir, "downloaded_matrices.csv"))
  failures <- read_csv_if_exists(file.path(out_dir, "failures.csv"))
  rows <- list()
  if (!is.null(success) && nrow(success) && "SeriesAccession" %in% names(success)) {
    rows[[length(rows) + 1L]] <- success
  }
  if (!is.null(failures) && nrow(failures) && "SeriesAccession" %in% names(failures)) {
    rows[[length(rows) + 1L]] <- failures
  }
  if (!length(rows)) {
    return(empty_bronze_manifest())
  }
  combined <- dplyr::bind_rows(rows)
  if (!"Status" %in% names(combined)) combined$Status <- NA_character_
  if (!"MatrixFile" %in% names(combined)) combined$MatrixFile <- NA_character_
  combined
}

summarize_pipeline_dir <- function(out_dir, n_datasets = NA_integer_, n_samples = NA_integer_) {
  success <- read_csv_if_exists(file.path(out_dir, "downloaded_matrices.csv"))
  failures <- read_csv_if_exists(file.path(out_dir, "failures.csv"))
  summary_csv <- read_csv_if_exists(file.path(out_dir, "summary.csv"))
  n_ok <- if (is.null(success)) 0L else nrow(success)
  n_fail <- if (is.null(failures)) 0L else nrow(failures)
  gse_ok <- if (n_ok && "SeriesAccession" %in% names(success)) {
    length(unique(success$SeriesAccession))
  } else {
    0L
  }
  out <- data.frame(
    Metric = c(
      "Downloaded matrix rows",
      "Unique downloaded studies",
      "Failure rows",
      "Datasets in DuckDB",
      "Samples in DuckDB"
    ),
    Value = c(n_ok, gse_ok, n_fail, n_datasets, n_samples),
    stringsAsFactors = FALSE
  )
  if (!is.null(summary_csv) && all(c("Metric", "Value") %in% names(summary_csv))) {
    out <- rbind(out, summary_csv[, c("Metric", "Value")])
  }
  out
}

write_filtered_metadata_xlsx <- function(plan, path, max_n = Inf) {
  if (!nrow(plan)) {
    stop("No studies match the current filters.")
  }
  keep <- head(seq_len(nrow(plan)), max_n)
  metadata <- data.frame(
    SeriesAccession = plan$GSE[keep],
    SeriesTitle = plan$StudyTitle[keep],
    Organism = plan$Organism[keep],
    stringsAsFactors = FALSE
  )
  writexl::write_xlsx(list(Metadata = metadata), path)
  invisible(path)
}

rscript_bin <- function() {
  found <- Sys.which("Rscript")
  if (nzchar(found)) {
    return(unname(found))
  }
  win <- "C:/Program Files/R/R-4.5.2/bin/Rscript.exe"
  if (file.exists(win)) {
    return(win)
  }
  "Rscript"
}

python_bin <- function() {
  found <- Sys.which("python")
  if (nzchar(found)) unname(found) else "python"
}
