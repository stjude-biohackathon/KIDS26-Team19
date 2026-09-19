#!/usr/bin/env Rscript
# Gold layer: harmonise treatment labels, score samples, rank agents.
#
# This is the step AGENTS.md describes as the point of the project: "scores samples
# with LSC6 and LSC17 and ranks drugs for lab follow-up". It reads series matrices,
# pulls a drug/vehicle label out of the free-text sample characteristics, scores every
# sample with pLSC6, and ranks agents by how far they move the stemness score away
# from their own study's control arm.
#
#   Rscript R/rank_drugs.R [dir_of_series_matrices] [signature]
#
# Cross-study comparability: the score is a weighted sum of expression, so raw values
# are not comparable across platforms or normalisations. Every contrast is therefore
# computed in within-study standard deviations, never in raw score units.

suppressPackageStartupMessages({ library(survival) })
source(file.path("R", "lsc_scores.R"))

args   <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) >= 1) args[1] else "Toy-Datasets"
SIG    <- if (length(args) >= 2) args[2] else "pLSC6"
# Output prefix. Defaults to the pilot names; pass a prefix to keep corpus-scale runs
# from overwriting them (e.g. "corpus_" -> results/corpus_drug_ranking.csv).
PFX    <- if (length(args) >= 3) args[3] else ""
OUT    <- "results"
on <- function(name) file.path(OUT, paste0(PFX, name))

source(file.path("R", "harmonise.R"))

# ---- per-study scoring -----------------------------------------------------

score_study <- function(path) {
  sm <- read_series_matrix(path)
  gse <- sm$gse

  # Both published scores assume log-scale input. Transform when the matrix is
  # plainly linear, and record that it happened.
  transformed <- FALSE
  if (isTRUE(max(sm$expr, na.rm = TRUE) > 50)) {
    sm$expr <- log2(pmax(sm$expr, 0) + 1)
    transformed <- TRUE
  }

  res <- try(suppressWarnings(score_lsc(sm, SIG)), silent = TRUE)
  if (inherits(res, "try-error")) {
    message(sprintf("  %-10s skipped: %s", gse,
                    substr(sub(".*: ", "", attr(res, "condition")$message), 1, 60)))
    return(NULL)
  }

  h <- harmonise(res$Characteristics)
  out <- cbind(res[, c("GSM", "SampleTitle", "Score", "GSE", "Platform")], h)

  if (all(is.na(out$CellLine))) {
    hit <- regmatches(out$SampleTitle,
                      regexpr(paste(KNOWN_LINES, collapse = "|"), out$SampleTitle,
                              ignore.case = TRUE))
    if (length(hit) == nrow(out) && all(nzchar(hit))) out$CellLine <- toupper(hit)
  }
  out$Log2Applied <- transformed

  if (all(is.na(out$Agent))) {
    message(sprintf("  %-10s scored (n=%d) but has no treatment field", gse, nrow(out)))
    return(out)
  }
  message(sprintf("  %-10s scored (n=%d), arms: %s", gse, nrow(out),
                  paste(sort(unique(na.omit(out$Agent))), collapse = ", ")))
  out
}

# ---- run -------------------------------------------------------------------

files <- list.files(in_dir, pattern = "series_matrix\\.txt\\.gz$", full.names = TRUE)
if (!length(files)) stop("No series matrices under ", in_dir, call. = FALSE)
message(sprintf("Scoring %d studies with %s\n", length(files), SIG))

scored <- do.call(rbind, lapply(files, score_study))
if (is.null(scored)) stop("Nothing scored.", call. = FALSE)

# Within-study z, so contrasts are comparable across platforms and normalisations.
scored$ScoreZ <- ave(scored$Score, scored$GSE, FUN = function(v) {
  s <- stats::sd(v, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(NA_real_, length(v)) else as.numeric(scale(v))
})

dir.create(OUT, showWarnings = FALSE)
write.csv(scored, on("sample_scores.csv"), row.names = FALSE)

# ---- rank ------------------------------------------------------------------

# A contrast needs a control arm in the same study, and a cell line where one is
# recorded, so treated and control cells are otherwise alike.
usable <- scored[!is.na(scored$Agent), , drop = FALSE]
usable$Strat <- ifelse(!is.na(usable$CellLine), usable$CellLine,
                ifelse(!is.na(usable$Genotype), usable$Genotype, NA_character_))
usable$Stratum <- ifelse(is.na(usable$Strat), usable$GSE,
                         paste(usable$GSE, usable$Strat, sep = " / "))

rank_rows <- list()
for (st in unique(usable$Stratum)) {
  d <- usable[usable$Stratum == st, ]
  ctrl <- d[d$Agent == "control", ]
  if (!nrow(ctrl)) next
  for (ag in setdiff(unique(d$Agent), "control")) {
    tr <- d[d$Agent == ag, ]
    if (!nrow(tr)) next
    rank_rows[[length(rank_rows) + 1]] <- data.frame(
      Stratum = st, GSE = d$GSE[1], CellLine = d$Strat[1], Agent = ag,
      Context = { g <- unique(na.omit(d$Genotype))
                  if (length(g)) g[1]
                  else if (grepl("NfsB|nitroreduct", d$Strat[1], ignore.case = TRUE)) "NfsB-positive"
                  else if (grepl("wild.?type", d$Strat[1], ignore.case = TRUE)) "wild type"
                  else NA_character_ },
      nTreated = nrow(tr), nControl = nrow(ctrl),
      DeltaZ = mean(tr$ScoreZ, na.rm = TRUE) - mean(ctrl$ScoreZ, na.rm = TRUE),
      Confidence = tr$AgentConfidence[1], Log2Applied = tr$Log2Applied[1],
      stringsAsFactors = FALSE)
  }
}

if (!length(rank_rows)) stop("No study had both a treated and a control arm.", call. = FALSE)
rk <- do.call(rbind, rank_rows)
rk <- rk[order(rk$DeltaZ), ]
rk$DeltaZ <- round(rk$DeltaZ, 3)
write.csv(rk, on("drug_ranking.csv"), row.names = FALSE)

cat("\n=== Agents ranked by shift in", SIG, "versus their own study control ===\n")
cat("    negative = drug lowers the stemness score (the direction of interest)\n\n")
print(rk[, c("GSE", "CellLine", "Agent", "nTreated", "nControl", "DeltaZ", "Confidence")],
      row.names = FALSE)

agg <- aggregate(DeltaZ ~ Agent, data = rk, FUN = mean)
agg$nContrasts <- as.integer(table(rk$Agent)[agg$Agent])
agg <- agg[order(agg$DeltaZ), ]
cat("\n=== Collapsed across strata ===\n\n")
print(agg, row.names = FALSE)

png(on("drug_ranking.png"), width = 1500, height = 950, res = 175)
op <- par(mar = c(5, 13, 4, 2))
cols <- ifelse(agg$DeltaZ < 0, "#2D6A4F", "#A8202B")
bp <- barplot(agg$DeltaZ, horiz = TRUE, names.arg = agg$Agent, las = 1, col = cols,
              border = NA, xlab = sprintf("mean shift in %s (within-study SD)", SIG),
              main = sprintf("Agents ranked by %s shift vs control", SIG),
              cex.names = 0.85)
abline(v = 0, col = "#5A6475", lwd = 1.4)
mtext(sprintf("%d contrasts across %d studies; negative favours loss of stemness",
              nrow(rk), length(unique(rk$GSE))), side = 3, line = 0.2, cex = 0.78,
      col = "#5A6475")
par(op); dev.off()

cat(sprintf("\nWrote %s, %s, %s\n",
            on("sample_scores.csv"), on("drug_ranking.csv"), on("drug_ranking.png")))
