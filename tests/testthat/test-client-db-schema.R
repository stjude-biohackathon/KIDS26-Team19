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

  run("initDb.py", c("--db-path", db_file, "--diagnosis", "aml"))
  run("updateDb.py", c(matrices_dir, "--diagnosis", "aml", "--db-path", db_file))

  list(work = work, db_file = db_file)
}

test_that("initDb and updateDb create sample table with treatment and control", {
  fixture <- build_temp_database()
  on.exit(unlink(fixture$work, recursive = TRUE), add = TRUE)

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = fixture$db_file, read_only = TRUE)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE, after = FALSE)

  tables <- DBI::dbListTables(con)
  expect_true(all(c("diagnosis", "dataset", "sample") %in% tables))
  expect_false("samples" %in% tables)

  sample_cols <- DBI::dbGetQuery(con, "DESCRIBE sample")$column_name
  expect_true("treatment" %in% sample_cols)
  expect_true("control" %in% sample_cols)
  expect_false("sample_characteristics_ch1" %in% sample_cols)
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

test_that("updateDb populates treatment and control from characteristics row", {
  py <- find_python_bin()
  if (!nzchar(py)) {
    testthat::skip("python3 not available")
  }
  toy <- file.path(repo_root, "Toy-Datasets", "GSE155640_series_matrix.txt.gz")
  txt <- file.path(repo_root, "Toy-Datasets", "GSE155640_series_matrix.txt")
  if (!file.exists(toy) && !file.exists(txt)) {
    testthat::skip("GSE155640 matrix fixture missing")
  }

  work <- tempfile("geo_treatment_")
  dir.create(work)
  matrices_dir <- file.path(work, "matrices")
  dir.create(matrices_dir)
  if (file.exists(toy)) {
    file.copy(toy, file.path(matrices_dir, basename(toy)))
  } else {
    dest <- file.path(matrices_dir, "GSE155640_series_matrix.txt.gz")
    con_in <- file(txt, "r")
    con_out <- gzfile(dest, "w")
    writeLines(readLines(con_in, warn = FALSE), con_out)
    close(con_in)
    close(con_out)
  }
  db_file <- file.path(work, "treatment.db")
  system2(py, c(
    file.path(repo_root, "parsing", "initDb.py"),
    "--db-path", db_file, "--diagnosis", "aml"
  ))
  system2(py, c(
    file.path(repo_root, "parsing", "updateDb.py"),
    matrices_dir, "--diagnosis", "aml", "--db-path", db_file
  ))

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_file, read_only = TRUE)
  filled <- DBI::dbGetQuery(
    con,
    "SELECT COUNT(*) AS n FROM sample WHERE treatment IS NOT NULL"
  )$n
  expect_true(filled > 0L)
  dmso_treatment <- DBI::dbGetQuery(
    con,
    "SELECT COUNT(*) AS n FROM sample WHERE treatment LIKE '%DMSO%'"
  )$n
  expect_equal(dmso_treatment, 0L)
  dmso_control <- DBI::dbGetQuery(
    con,
    "SELECT COUNT(*) AS n FROM sample WHERE control LIKE '%DMSO%' OR control LIKE '%dmso%'"
  )$n
  expect_true(dmso_control > 0L)
  drug <- DBI::dbGetQuery(
    con,
    "SELECT COUNT(*) AS n FROM sample WHERE treatment LIKE '%OG86%'"
  )$n
  expect_true(drug > 0L)
  DBI::dbDisconnect(con, shutdown = TRUE)
  unlink(work, recursive = TRUE)
})
