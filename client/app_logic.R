## Data access for the Shiny app: reads the persistent DuckDB at data/geo.db.
## The database is created once by parsing/initDb.py and loaded by
## parsing/updateDb.py; this app only reads it.

default_db_path <- function(repo_root) file.path(repo_root, "data", "geo.db")

REQUIRED_TABLES <- c("diagnosis", "dataset", "sample")

SAMPLES_TABLE_COLUMNS <- c(
  "diagnosis_name",
  "series_accession",
  "series_platform_id",
  "series_pubmed_id",
  "sample_geo_accession",
  "sample_organism_ch1",
  "sample_data_row_count",
  "treatment",
  "control",
  "drug",
  "sample_molecule_ch1"
)

samples_table_column_labels <- c(
  "Diagnosis",
  "Series",
  "Platform",
  "PubMed ID",
  "Sample",
  "Organism",
  "Row count",
  "Treatment",
  "Control",
  "Drug",
  "Sample molecule"
)

SAMPLES_TABLE_SELECT <- paste(
  "SELECT",
  "g.diagnosis_name,",
  "d.series_geo_accession AS series_accession,",
  "d.series_platform_id,",
  "d.series_pubmed_id,",
  "s.sample_geo_accession,",
  "s.sample_organism_ch1,",
  "s.sample_data_row_count,",
  "s.treatment,",
  "s.control,",
  "s.drug,",
  "s.sample_molecule_ch1",
  "FROM sample s",
  "JOIN dataset d ON d.dataset_id = s.dataset_id",
  "JOIN diagnosis g ON g.diagnosis_id = s.diagnosis_id"
)

SAMPLES_TABLE_ORDER <-
  "ORDER BY diagnosis_name, series_geo_accession, sample_geo_accession"

samples_table_schema_ok <- function(connection) {
  cols <- DBI::dbGetQuery(connection, "DESCRIBE sample")$column_name
  required <- c(
    "sample_geo_accession",
    "treatment",
    "control",
    "drug",
    "sample_molecule_ch1"
  )
  missing <- setdiff(required, cols)
  list(ok = length(missing) == 0L, columns = cols, missing = missing)
}

#' Open the persistent database read-only, or NULL if it is missing or unbuilt.
connect_geo_database <- function(db_path) {
  if (!file.exists(db_path)) {
    return(NULL)
  }
  con <- tryCatch(
    DBI::dbConnect(duckdb::duckdb(), dbdir = db_path, read_only = TRUE),
    error = function(e) NULL
  )
  if (is.null(con)) {
    return(NULL)
  }
  ok <- tryCatch(
    all(REQUIRED_TABLES %in% DBI::dbListTables(con)),
    error = function(e) FALSE
  )
  if (!isTRUE(ok)) {
    DBI::dbDisconnect(con, shutdown = TRUE)
    return(NULL)
  }
  schema <- samples_table_schema_ok(con)
  if (!schema$ok) {
    DBI::dbDisconnect(con, shutdown = TRUE)
    return(NULL)
  }
  con
}

#' Diagnosis names present in the database, in display order.
list_diagnoses <- function(connection) {
  DBI::dbGetQuery(
    connection,
    "SELECT diagnosis_name FROM diagnosis ORDER BY diagnosis_name"
  )$diagnosis_name
}

#' Sample rows joined to dataset and diagnosis; diagnosis = NULL returns every row.
fetch_samples_table <- function(connection, diagnosis = NULL) {
  if (is.null(diagnosis) || !nzchar(diagnosis)) {
    return(DBI::dbGetQuery(
      connection,
      paste(SAMPLES_TABLE_SELECT, SAMPLES_TABLE_ORDER)
    ))
  }
  DBI::dbGetQuery(
    connection,
    paste(SAMPLES_TABLE_SELECT, "WHERE g.diagnosis_name = ?", SAMPLES_TABLE_ORDER),
    params = list(tolower(diagnosis))
  )
}

missing_database_message <- function(db_path) {
  paste0(
    "No usable database at ", db_path, ".\n",
    "Create it once with:  python parsing/initDb.py --diagnosis aml\n",
    "Then load matrices:   python parsing/updateDb.py <matrices dir> --diagnosis aml\n",
    "The database is committed to the repo, so `git pull` usually supplies it."
  )
}
