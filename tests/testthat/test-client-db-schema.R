source(file.path("..", "..", "client", "app_logic.R"))

repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))

find_python_bin <- function() {
  venv <- file.path(repo_root, "parsing", ".venv", "bin", "python3")
  if (file.exists(venv)) venv else Sys.which("python3")
}

build_temp_database <- function() {
  py <- find_python_bin()
  if (!nzchar(py)) {
    testthat::skip("python3 not available")
  }
  toy_matrix <- file.path(repo_root, "Toy-Datasets", "GSE128103_series_matrix.txt.gz")
  if (!file.exists(toy_matrix)) {
    testthat::skip("Toy-Datasets matrix fixture missing")
  }

  work <- tempfile("geo_db_test_")
  dir.create(work)
  matrices_dir <- file.path(work, "matrices")
  dir.create(matrices_dir)
  file.copy(toy_matrix, file.path(matrices_dir, basename(toy_matrix)))
  db_file <- file.path(work, "geo.db")

  run <- function(script, args) {
    out <- suppressWarnings(system2(
      py, c(file.path(repo_root, "parsing", script), args), stdout = TRUE, stderr = TRUE
    ))
    code <- attr(out, "status")
    if (!is.null(code) && code != 0L) {
      testthat::fail(paste(c(script, out), collapse = "\n"))
    }
    out
  }

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_file, read_only = TRUE)
  sample_cols <- DBI::dbGetQuery(con, "DESCRIBE samples")$column_name
  expect_true("series_accession" %in% sample_cols)
  expect_true("series_pubmed_id" %in% sample_cols)
  expect_true("sample_geo_accession" %in% sample_cols)
  expect_true("sample_characteristics_ch1" %in% sample_cols)
  expect_true("sample_molecule_ch1" %in% sample_cols)
  expect_false("SeriesAccession" %in% sample_cols)

  list(work = work, db_file = db_file)
}

test_that("initDb creates the diagnosis, dataset and sample tables", {
  fixture <- build_temp_database()
  on.exit(unlink(fixture$work, recursive = TRUE), add = TRUE)

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = fixture$db_file, read_only = TRUE)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE, after = FALSE)

  tables <- DBI::dbListTables(con)
  expect_true(all(c("diagnosis", "dataset", "sample") %in% tables))
  expect_false("samples" %in% tables)

  dataset_cols <- DBI::dbGetQuery(con, "DESCRIBE dataset")$column_name
  expect_true(all(
    c("dataset_id", "diagnosis_id", "source_file", "series_geo_accession") %in% dataset_cols
  ))

  sample_cols <- DBI::dbGetQuery(con, "DESCRIBE sample")$column_name
  expect_true(all(
    c("sample_id", "dataset_id", "diagnosis_id", "sample_geo_accession") %in% sample_cols
  ))
})

test_that("connect_geo_database opens the persistent database and lists diagnoses", {
  fixture <- build_temp_database()
  on.exit(unlink(fixture$work, recursive = TRUE), add = TRUE)

  con <- connect_geo_database(fixture$db_file)
  expect_false(is.null(con))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE, after = FALSE)

  expect_identical(list_diagnoses(con), "aml")
})

test_that("connect_geo_database returns NULL for a missing or unrelated database", {
  expect_null(connect_geo_database(file.path(tempdir(), "does-not-exist.db")))

  stray <- tempfile("stray_", fileext = ".db")
  on.exit(unlink(stray), add = TRUE)
  stray_con <- DBI::dbConnect(duckdb::duckdb(), dbdir = stray)
  DBI::dbExecute(stray_con, "CREATE TABLE unrelated (x INTEGER)")
  DBI::dbDisconnect(stray_con, shutdown = TRUE)

  expect_null(connect_geo_database(stray))
})
