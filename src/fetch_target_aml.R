#!/usr/bin/env Rscript
#
# Pull open-access TARGET-AML clinical data from the NCI GDC and produce a
# Kaplan-Meier baseline.
#
# Why this cohort: TARGET-AML exposes vital status, days to death, and days to
# last follow-up through the public GDC API with no authentication and no data
# access request. That gives patient-level survival labels on the same subjects
# that have open gene expression, somatic mutation, and methylation files, which
# is what an outcome-prediction model needs. It is also a pediatric cohort.
#
# Nothing here downloads controlled-access data. Everything written to disk is
# derived from open files and is safe to commit.
#
# Usage:
#   Rscript src/fetch_target_aml.R                  # default TARGET-AML
#   Rscript src/fetch_target_aml.R TCGA-LAML        # any GDC project id
#
# Requires: jsonlite, survival. Optional: ggplot2 (falls back to base graphics).

suppressPackageStartupMessages({
  library(jsonlite)
  library(survival)
})

GDC_CASES   <- "https://api.gdc.cancer.gov/cases"
PAGE_SIZE   <- 500L
CACHE_DIR   <- "data/processed"
RESULTS_DIR <- "results"

CLINICAL_FIELDS <- paste(
  "submitter_id",
  "demographic.vital_status",
  "demographic.days_to_death",
  "demographic.gender",
  "demographic.race",
  "diagnoses.age_at_diagnosis",
  "diagnoses.days_to_last_follow_up",
  "diagnoses.primary_diagnosis",
  sep = ","
)

# ---- helpers ---------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x

# Pull a single page of cases. GDC omits null fields entirely, so every
# accessor below has to tolerate an absent element rather than assume a shape.
gdc_page <- function(project_id, from, size) {
  body <- list(
    filters = list(
      op = "in",
      content = list(field = "project.project_id", value = list(project_id))
    ),
    fields = CLINICAL_FIELDS,
    format = "json",
    size   = as.character(size),
    from   = as.character(from)
  )
  res <- tryCatch(
    fromJSON(
      rawToChar(curl_post(GDC_CASES, toJSON(body, auto_unbox = TRUE))),
      simplifyVector = FALSE
    ),
    error = function(e) stop("GDC request failed: ", conditionMessage(e), call. = FALSE)
  )
  res$data
}

# Minimal POST without pulling in httr, so the script runs on a bare R install.
curl_post <- function(url, json) {
  tmp <- tempfile(fileext = ".json")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(json, tmp)
  out <- system2(
    "curl",
    c("-s", "--max-time", "120", "-H", shQuote("Content-Type: application/json"),
      "-X", "POST", shQuote(url), "--data-binary", paste0("@", tmp)),
    stdout = TRUE
  )
  charToRaw(paste(out, collapse = ""))
}

fetch_clinical <- function(project_id) {
  first <- gdc_page(project_id, from = 0L, size = 1L)
  total <- first$pagination$total
  if (is.null(total) || total == 0L) {
    stop("No cases returned for project '", project_id, "'.", call. = FALSE)
  }
  message(sprintf("%s: %d cases reported by GDC", project_id, total))

  hits <- list()
  from <- 0L
  while (from < total) {
    page <- gdc_page(project_id, from = from, size = PAGE_SIZE)
    hits <- c(hits, page$hits)
    from <- from + PAGE_SIZE
    message(sprintf("  fetched %d / %d", min(from, total), total))
  }
  hits
}

# Flatten one case record. Takes the first diagnosis when several are present;
# TARGET-AML records a small number of subjects with more than one.
flatten_case <- function(h) {
  dg <- if (length(h$diagnoses)) h$diagnoses[[1]] else list()
  dm <- h$demographic %||% list()
  data.frame(
    submitter_id          = h$submitter_id            %||% NA_character_,
    vital_status          = dm$vital_status           %||% NA_character_,
    days_to_death         = as.numeric(dm$days_to_death %||% NA),
    gender                = dm$gender                 %||% NA_character_,
    race                  = dm$race                   %||% NA_character_,
    age_at_diagnosis_days = as.numeric(dg$age_at_diagnosis %||% NA),
    days_to_last_followup = as.numeric(dg$days_to_last_follow_up %||% NA),
    primary_diagnosis     = dg$primary_diagnosis      %||% NA_character_,
    stringsAsFactors      = FALSE
  )
}

# Standard overall-survival encoding: time runs to death when the subject died,
# otherwise to last known follow-up, and the event flag marks death.
build_survival_table <- function(df) {
  df$time  <- ifelse(df$vital_status == "Dead", df$days_to_death, df$days_to_last_followup)
  df$event <- as.integer(df$vital_status == "Dead")
  df$age_years <- df$age_at_diagnosis_days / 365.25

  # Age bands follow how pediatric AML risk is usually discussed: infants do
  # worse, and outcomes decline again in adolescents and young adults.
  df$age_group <- cut(
    df$age_years,
    breaks = c(-Inf, 1, 10, 15, Inf),
    labels = c("Infant (<1y)", "Child (1-9y)", "Adolescent (10-14y)", "AYA (15y+)"),
    right  = FALSE
  )

  keep <- !is.na(df$time) & !is.na(df$event) & df$time >= 0
  dropped <- sum(!keep)
  if (dropped > 0) {
    message(sprintf("  dropped %d of %d cases with no usable follow-up time",
                    dropped, nrow(df)))
  }
  df[keep, , drop = FALSE]
}

plot_km <- function(surv_df, project_id, outfile) {
  fit <- survfit(Surv(time / 365.25, event) ~ age_group, data = surv_df)
  png(outfile, width = 1600, height = 1100, res = 180)
  on.exit(dev.off(), add = TRUE)
  plot(
    fit,
    col  = c("#2B6CB0", "#2F855A", "#B7791F", "#9B2C2C"),
    lwd  = 2.2,
    xlab = "Years from diagnosis",
    ylab = "Overall survival",
    mark.time = TRUE,
    main = sprintf("%s overall survival by age at diagnosis (n = %d)",
                   project_id, nrow(surv_df))
  )
  legend(
    "bottomleft",
    legend = levels(surv_df$age_group),
    col    = c("#2B6CB0", "#2F855A", "#B7791F", "#9B2C2C"),
    lwd    = 2.2, bty = "n", cex = 0.85
  )
  invisible(fit)
}

# ---- main ------------------------------------------------------------------

main <- function(project_id = "TARGET-AML") {
  dir.create(CACHE_DIR,   recursive = TRUE, showWarnings = FALSE)
  dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

  slug        <- tolower(gsub("[^A-Za-z0-9]+", "_", project_id))
  clinical_fp <- file.path(CACHE_DIR,   paste0(slug, "_clinical.csv"))
  km_fp       <- file.path(RESULTS_DIR, paste0(slug, "_km_by_age.png"))

  # Re-runnable: a second run reuses the cached pull instead of hitting the API.
  if (file.exists(clinical_fp)) {
    message("Reusing cached clinical table: ", clinical_fp)
    clin <- read.csv(clinical_fp, stringsAsFactors = FALSE)
  } else {
    hits <- fetch_clinical(project_id)
    clin <- do.call(rbind, lapply(hits, flatten_case))
    write.csv(clin, clinical_fp, row.names = FALSE)
    message("Wrote ", clinical_fp)
  }

  surv_df <- build_survival_table(clin)

  message("\n--- cohort ---")
  message(sprintf("  cases with usable follow-up : %d", nrow(surv_df)))
  message(sprintf("  deaths                      : %d (%.1f%%)",
                  sum(surv_df$event), 100 * mean(surv_df$event)))
  message(sprintf("  median follow-up (years)    : %.2f",
                  median(surv_df$time, na.rm = TRUE) / 365.25))
  message("  by age group:")
  print(table(surv_df$age_group, dnn = NULL))

  fit <- survfit(Surv(time / 365.25, event) ~ age_group, data = surv_df)
  message("\n--- 5-year overall survival ---")
  s5 <- summary(fit, times = 5, extend = TRUE)
  for (i in seq_along(s5$strata)) {
    message(sprintf("  %-22s %.3f (95%% CI %.3f-%.3f)",
                    sub("age_group=", "", as.character(s5$strata[i])),
                    s5$surv[i], s5$lower[i], s5$upper[i]))
  }

  lr <- survdiff(Surv(time, event) ~ age_group, data = surv_df)
  p  <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
  message(sprintf("\n  log-rank across age groups: chisq = %.1f, p = %.3g", lr$chisq, p))

  plot_km(surv_df, project_id, km_fp)
  message("\nWrote ", km_fp)

  write.csv(surv_df, file.path(CACHE_DIR, paste0(slug, "_survival.csv")),
            row.names = FALSE)
  message("Wrote ", file.path(CACHE_DIR, paste0(slug, "_survival.csv")))

  invisible(list(clinical = clin, survival = surv_df, fit = fit))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  main(if (length(args)) args[1] else "TARGET-AML")
}
