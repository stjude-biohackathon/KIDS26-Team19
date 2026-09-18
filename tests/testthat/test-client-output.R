source(file.path("..", "..", "client", "app_logic.R"))

test_that("samples query selects the joined columns and matches its labels", {
  expect_identical(length(SAMPLES_TABLE_COLUMNS), length(samples_table_column_labels))
  expect_true(grepl("FROM sample s", SAMPLES_TABLE_SELECT, fixed = TRUE))
  expect_true(grepl("JOIN dataset d", SAMPLES_TABLE_SELECT, fixed = TRUE))
  expect_true(grepl("JOIN diagnosis g", SAMPLES_TABLE_SELECT, fixed = TRUE))
  expect_false(grepl("FROM samples", SAMPLES_TABLE_SELECT, fixed = TRUE))
})

test_that("default_db_path points at the persistent database", {
  expect_identical(default_db_path("/tmp/repo"), file.path("/tmp/repo", "data", "geo.db"))
})

test_that("missing_database_message names the init and update scripts", {
  message_text <- missing_database_message("data/geo.db")
  expect_match(message_text, "parsing/initDb.py", fixed = TRUE)
  expect_match(message_text, "parsing/updateDb.py", fixed = TRUE)
})