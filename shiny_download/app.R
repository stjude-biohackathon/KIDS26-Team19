## GEO Download and Sample Explorer
## Run from the kp_d3 repo root:  shiny::runApp("shiny_download")

library(shiny)
library(bslib)
library(DT)
library(DBI)
library(duckdb)

repo_root <- normalizePath(
  if (file.exists("parsing/initDb.py")) "." else "..",
  winslash = "/",
  mustWork = TRUE
)
planner_file <- file.path(repo_root, "shiny_download", "planner_logic.R")
if (!file.exists(planner_file)) planner_file <- "planner_logic.R"
source(planner_file, local = FALSE)

script05 <- file.path(repo_root, "R_Scripts", "05_Download_Metadata_Inventory.R")
default_workbook <- default_heme_inventory_path(repo_root)
default_out_dir <- file.path(repo_root, "downloads", "shiny_download")
default_db_path <- file.path(repo_root, "data", "geo.db")

all_choice <- c("All" = "")

dropdown <- function(id, label) {
  selectInput(id, label, choices = all_choice, selected = "", selectize = FALSE, width = "100%")
}

ui <- page_sidebar(
  title = "GEO Download and Sample Explorer",
  theme = bs_theme(version = 5, bootswatch = "flatly", font_scale = 0.95),
  sidebar = sidebar(
    width = 340,
    h5("Inventory"),
    textInput("workbook", "Excel workbook", default_workbook),
    fileInput("inventory_file", "Upload Excel", accept = c(".xlsx", ".xls")),
    actionButton("reload_inventory", "Load inventory", class = "btn-primary w-100"),
    hr(),
    h5("Filters"),
    dropdown("diseases", "Disease"),
    dropdown("cell_lines", "Cell line"),
    dropdown("organisms", "Organism"),
    dropdown("specimen_types", "Specimen type"),
    dropdown("treatment_statuses", "Treatment status"),
    checkboxInput("require_vehicle", "Only DMSO / vehicle-control series", FALSE),
    hr(),
    h5("Run"),
    textInput("query", "Keyword query (optional)", ""),
    textInput("out_dir", "Output directory", default_out_dir),
    textInput("db_path", "DuckDB path", default_db_path),
    textInput("diagnosis", "Diagnosis name", "heme_cancer"),
    numericInput("max_studies", "Max studies to download", 5, min = 1, step = 1),
    actionButton("run", "Download, update DB, and explore", class = "btn-warning w-100")
  ),
  layout_columns(
    col_widths = c(3, 3, 3, 3),
    fill = FALSE,
    value_box("Matching studies", textOutput("n_studies"), theme = "primary"),
    value_box("Matching samples", textOutput("n_samples"), theme = "primary"),
    value_box("Downloaded files", textOutput("n_downloaded"), theme = "success"),
    value_box("DuckDB samples", textOutput("n_db_samples"), theme = "secondary")
  ),
  p(class = "text-muted mb-2", textOutput("status", inline = TRUE)),
  navset_card_underline(
    id = "main_tabs",
    nav_panel("Matching studies", DTOutput("plan_table")),
    nav_panel("Matching samples", DTOutput("sample_table")),
    nav_panel("Downloaded", DTOutput("downloaded")),
    nav_panel("Failures", DTOutput("failures")),
    nav_panel("Datasets", DTOutput("datasets")),
    nav_panel("Samples", DTOutput("samples")),
    nav_panel("Run log", verbatimTextOutput("log"))
  )
)

server <- function(input, output, session) {
  status <- reactiveVal("Ready.")
  logtxt <- reactiveVal("")
  version <- reactiveVal(0L)
  con <- reactiveVal(NULL)
  inventory_data <- reactiveVal(NULL)

  refresh_db <- function() {
    old <- isolate(con())
    if (!is.null(old)) try(dbDisconnect(old, shutdown = TRUE), silent = TRUE)
    path <- isolate(input$db_path)
    con(if (nzchar(path) && file.exists(path)) {
      dbConnect(duckdb(shared_home = FALSE), path, read_only = TRUE)
    } else {
      NULL
    })
    version(isolate(version()) + 1L)
  }

  observeEvent(input$db_path, refresh_db(), ignoreNULL = TRUE)

  set_dropdown <- function(id, values) {
    updateSelectInput(
      session,
      id,
      choices = c(all_choice, stats::setNames(values, values)),
      selected = ""
    )
  }

  load_inventory <- function(path) {
    tryCatch(
      {
        data <- read_heme_inventory(path)
        inventory_data(data)
        choices <- heme_filter_choices(data)
        set_dropdown("diseases", choices$diseases)
        set_dropdown("cell_lines", choices$cell_lines)
        set_dropdown("organisms", choices$organisms)
        set_dropdown("specimen_types", choices$specimen_types)
        set_dropdown("treatment_statuses", choices$treatment_statuses)
        status(sprintf("Loaded %s samples / %s series from %s",
                       nrow(data), length(unique(data$GSE)), basename(path)))
      },
      error = function(e) {
        inventory_data(NULL)
        status(paste("Could not read inventory:", conditionMessage(e)))
      }
    )
  }

  if (file.exists(default_workbook)) load_inventory(default_workbook)

  observeEvent(input$reload_inventory, {
    path <- input$workbook
    if (!is.null(input$inventory_file) && nzchar(input$inventory_file$datapath)) {
      path <- input$inventory_file$datapath
    }
    load_inventory(path)
  })

  selected_or_null <- function(value) {
    if (is.null(value) || !nzchar(value)) NULL else value
  }

  filtered_samples <- reactive({
    data <- inventory_data()
    if (is.null(data)) return(data.frame())
    filter_heme_inventory(
      data,
      diseases = selected_or_null(input$diseases),
      cell_lines = selected_or_null(input$cell_lines),
      organisms = selected_or_null(input$organisms),
      specimen_types = selected_or_null(input$specimen_types),
      treatment_statuses = selected_or_null(input$treatment_statuses),
      require_vehicle_control = isTRUE(input$require_vehicle)
    )
  })

  pipeline_tables <- reactive({
    version()
    dir <- input$out_dir
    list(
      downloaded = read_csv_if_exists(file.path(dir, "downloaded_matrices.csv")),
      failures = read_csv_if_exists(file.path(dir, "failures.csv"))
    )
  })

  current_plan <- reactive({
    samples <- filtered_samples()
    if (!nrow(samples)) return(heme_download_plan(samples))
    manifest <- pipeline_tables()$downloaded
    if (is.null(manifest)) manifest <- empty_bronze_manifest()
    merge_plan_with_downloads(heme_download_plan(samples), manifest)
  })

  db_n <- reactive({
    version()
    x <- con()
    if (is.null(x)) return(list(datasets = 0L, samples = 0L))
    list(
      datasets = tryCatch(dbGetQuery(x, "SELECT COUNT(*) AS n FROM dataset")$n, error = function(e) 0L),
      samples = tryCatch(dbGetQuery(x, "SELECT COUNT(*) AS n FROM sample")$n, error = function(e) 0L)
    )
  })

  output$status <- renderText(status())
  output$n_studies <- renderText({
    plan <- current_plan()
    if (!nrow(plan)) "0" else as.character(nrow(plan))
  })
  output$n_samples <- renderText(as.character(nrow(filtered_samples())))
  output$n_downloaded <- renderText({
    d <- pipeline_tables()$downloaded
    as.character(if (is.null(d)) 0L else nrow(d))
  })
  output$n_db_samples <- renderText(as.character(db_n()$samples))

  show_table <- function(tbl) {
    datatable(
      tbl,
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 15, scrollX = TRUE, autoWidth = TRUE)
    )
  }

  output$plan_table <- renderDT({
    validate(need(nrow(current_plan()) > 0, "Load an inventory to see matching studies."))
    keep <- intersect(
      c("GSE", "Diseases", "StudyTitle", "Organism", "CellLines",
        "MatchingSamples", "MatrixFileCount", "AlreadyDownloaded", "DownloadStatus"),
      names(current_plan())
    )
    show_table(current_plan()[, keep, drop = FALSE])
  })

  output$sample_table <- renderDT({
    tbl <- filtered_samples()
    validate(need(nrow(tbl) > 0, "Load an inventory to see matching samples."))
    keep <- intersect(
      c("CancerType", "GSE", "GSM", "SampleTitle", "Organism", "SpecimenType",
        "TreatmentStatus", "TreatmentAgents", "CellLineNames", "Platform"),
      names(tbl)
    )
    show_table(tbl[, keep, drop = FALSE])
  })

  output$downloaded <- renderDT({
    tbl <- pipeline_tables()$downloaded
    validate(need(!is.null(tbl) && nrow(tbl) > 0, "No downloaded matrices yet."))
    show_table(tbl)
  })
  output$failures <- renderDT({
    tbl <- pipeline_tables()$failures
    validate(need(!is.null(tbl) && nrow(tbl) > 0, "No failures recorded."))
    show_table(tbl)
  })

  db_query <- function(sql) {
    x <- con()
    validate(need(!is.null(x), "No DuckDB file found. Run a download or check the DuckDB path."))
    dbGetQuery(x, sql)
  }
  output$datasets <- renderDT({
    version()
    show_table(db_query("SELECT * FROM dataset ORDER BY dataset_id"))
  })
  output$samples <- renderDT({
    version()
    show_table(db_query("SELECT * FROM sample ORDER BY sample_id"))
  })
  output$log <- renderText({
    txt <- logtxt()
    if (!nzchar(txt)) "No run log yet." else txt
  })

  observeEvent(input$run, {
    useq <- nzchar(trimws(input$query))
    dir.create(input$out_dir, recursive = TRUE, showWarnings = FALSE)
    workbook <- input$workbook
    if (!is.null(input$inventory_file) && nzchar(input$inventory_file$datapath)) {
      workbook <- input$inventory_file$datapath
    }

    if (!useq) {
      plan <- current_plan()
      if (nrow(plan)) {
        workbook <- file.path(input$out_dir, "filtered_metadata.xlsx")
        write_filtered_metadata_xlsx(plan, workbook, input$max_studies)
      } else if (!file.exists(workbook)) {
        status("Workbook not found.")
        return()
      }
    }

    status("Running GEO download...")
    args <- if (useq) {
      c("--vanilla", script05, "--query", input$query, "--out", input$out_dir)
    } else {
      c("--vanilla", script05, workbook, input$out_dir)
    }
    r <- tryCatch(
      system2(rscript_bin(), args, stdout = TRUE, stderr = TRUE),
      error = function(e) paste("ERROR:", conditionMessage(e))
    )
    if (any(grepl("ERROR|Execution halted", r))) {
      logtxt(paste(r, collapse = "\n"))
      status("Download failed. See Run log.")
      return()
    }

    status("Updating database...")
    old <- isolate(con())
    if (!is.null(old)) {
      try(dbDisconnect(old, shutdown = TRUE), silent = TRUE)
      con(NULL)
    }
    if (!file.exists(input$db_path)) {
      system2(
        python_bin(),
        c(file.path(repo_root, "parsing", "initDb.py"),
          "--db-path", input$db_path, "--diagnosis", input$diagnosis),
        stdout = TRUE, stderr = TRUE
      )
    }
    u <- system2(
      python_bin(),
      c(file.path(repo_root, "parsing", "updateDb.py"),
        file.path(input$out_dir, "matrices"),
        "--db-path", input$db_path,
        "--diagnosis", input$diagnosis,
        "--report", file.path(input$out_dir, "db_update_report.tsv")),
      stdout = TRUE, stderr = TRUE
    )
    logtxt(paste(c(r, u), collapse = "\n"))
    status("Complete. Reports and database updated.")
    refresh_db()
  })

  session$onSessionEnded(function() {
    old <- tryCatch(isolate(con()), error = function(e) NULL)
    if (!is.null(old)) try(dbDisconnect(old, shutdown = TRUE), silent = TRUE)
  })
}

shinyApp(ui, server)
