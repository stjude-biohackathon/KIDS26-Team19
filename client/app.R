## Shiny app: browse GEO sample metadata in the persistent data/geo.db.
## R setup: install.packages(c("shiny", "DT", "DBI", "duckdb"))
## The database is created by parsing/initDb.py and loaded by parsing/updateDb.py;
## this app opens it read-only and never rebuilds it.

library(shiny)
library(DT)
library(DBI)
library(duckdb)

find_repo_root <- function() {
  for (candidate in c(".", "..")) {
    if (file.exists(file.path(candidate, "parsing", "initDb.py"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  stop("Cannot locate repository root (expected parsing/initDb.py).")
}

repo_root <- find_repo_root()
source(file.path(repo_root, "client", "app_logic.R"))
db_path <- default_db_path(repo_root)

ui <- fluidPage(
  titlePanel("GEO Sample Explorer"),
  sidebarLayout(
    sidebarPanel(
      textInput("workbook", "Inventory workbook", "/Users/fauzanahmed/Downloads/Heme_Cancer_GEO_Inventory.xlsx"),
      textInput("out_dir", "Download directory", file.path(repo_root, "downloads", "shiny_pipeline")),
      textInput("diagnosis_name", "Diagnosis name", "heme_cancer"),
      actionButton("run_pipeline", "Run download and add to database", class = "btn-primary"),
      tags$hr(),
      selectInput("diagnosis", "Diagnosis", choices = c("All diagnoses" = "")),
      helpText(
        "Read-only view of data/geo.db. Add studies outside this app with ",
        "parsing/updateDb.py. The database is preserved between runs."
      )
    ),
    mainPanel(
      h4(textOutput("state", inline = TRUE)),
      verbatimTextOutput("message"),
      DTOutput("results")
    )
  )
)

server <- function(input, output, session) {
  connection <- reactiveVal(connect_geo_database(db_path))
  status <- reactiveVal("")
  query_error <- reactiveVal("")

  if (is.null(connection())) {
    status(missing_database_message(db_path))
  } else {
    diagnoses <- list_diagnoses(connection())
    updateSelectInput(
      session,
      "diagnosis",
      choices = c("All diagnoses" = "", stats::setNames(diagnoses, toupper(diagnoses))),
      selected = ""
    )
    status(paste0("Opened ", db_path, " (", length(diagnoses), " diagnoses)."))
    session$onSessionEnded(function() dbDisconnect(connection(), shutdown = TRUE))
  }

  observeEvent(input$run_pipeline, {
    req(file.exists(input$workbook))
    dir.create(input$out_dir, recursive = TRUE, showWarnings = FALSE)
    output$message <- renderText("Running GEO download and database update...")
    r <- system2("Rscript", c(file.path(repo_root, "R_Scripts", "05_Download_Metadata_Inventory.R"), input$workbook, input$out_dir), stdout = TRUE, stderr = TRUE)
    system2("python", c(file.path(repo_root, "parsing", "updateDb.py"), file.path(input$out_dir, "matrices"), "--db-path", db_path, "--diagnosis", input$diagnosis_name, "--report", file.path(input$out_dir, "db_update_report.tsv")), stdout = TRUE, stderr = TRUE)
    output$message <- renderText(paste(c(r, "Database updated."), collapse = "\n"))
  })

  results <- reactive({
    req(!is.null(connection()))
    tryCatch(
      {
        rows <- fetch_samples_table(connection(), input$diagnosis)
        query_error("")
        rows
      },
      error = function(e) {
        query_error(conditionMessage(e))
        NULL
      }
    )
  })

  output$state <- renderText({
    if (is.null(connection())) {
      return("No database loaded.")
    }
    if (nzchar(query_error())) {
      return(paste0("Database query error: ", query_error()))
    }
    rows <- results()
    if (is.null(rows)) {
      return("Database loaded.")
    }
    paste0("Database loaded. ", nrow(rows), " sample rows.")
  })

  output$message <- renderText({ status() })

  output$results <- renderDT({
    validate(need(!is.null(connection()), missing_database_message(db_path)))
    tbl <- results()
    if (nzchar(query_error())) {
      validate(need(FALSE, paste0("Could not read samples:\n", query_error())))
    }
    validate(need(!is.null(tbl) && nrow(tbl) > 0L, "No sample rows to display."))
    treatment_idx <- which(names(tbl) == "treatment") - 1L
    column_defs <- list()
    if (length(treatment_idx) == 1L && treatment_idx >= 0L) {
      column_defs <- list(
        list(width = "320px", targets = treatment_idx)
      )
    }
    datatable(
      tbl,
      rownames = FALSE,
      filter = "top",
      colnames = samples_table_column_labels,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        columnDefs = column_defs
      )
    )
  })
}

shinyApp(ui, server)
