#!/usr/bin/env Rscript
# How many of the downloaded GEO series actually support the project's premise?
#
# The team's bronze layer holds 585 AML series. The pipeline needs studies with a
# treated arm and a vehicle or control arm. This reads the metadata header of every
# series, applies the same harmoniser the scoring code uses, and reports how many
# qualify. Header-only, so it costs no expression parsing.
#
#   Rscript R/scan_corpus.R

source(file.path("R", "harmonise.R"))

hdrs <- list.files(file.path("bronze", "headers"), pattern = "\\.hdr$", full.names = TRUE)
if (!length(hdrs)) stop("No headers under bronze/headers/.", call. = FALSE)
message(sprintf("Scanning %d series headers", length(hdrs)))

read_header <- function(fp) {
  ln <- readLines(fp, warn = FALSE)
  get <- function(key) {
    i <- grep(paste0("^!", key, "\t"), ln)
    if (!length(i)) return(NULL)
    lapply(i, function(j) gsub('^"|"$', "", strsplit(ln[j], "\t", fixed = TRUE)[[1]][-1]))
  }
  gsm <- get("Sample_geo_accession")
  if (is.null(gsm)) return(NULL)
  gsm <- gsm[[1]]
  ch <- get("Sample_characteristics_ch1")
  chars <- if (is.null(ch)) rep("", length(gsm)) else
    apply(do.call(rbind, lapply(ch, function(v) { length(v) <- length(gsm)
                                                  ifelse(is.na(v), "", v) })),
          2, function(col) paste(col[nzchar(col)], collapse = " | "))
  list(gse = (get("Series_geo_accession") %||% list(NA))[[1]][1],
       gpl = (get("Series_platform_id")   %||% list(NA))[[1]][1],
       title = (get("Series_title")       %||% list(NA))[[1]][1],
       org = (get("Sample_organism_ch1")  %||% list(NA))[[1]][1],
       n = length(gsm), chars = chars)
}

rows <- list()
for (fp in hdrs) {
  h <- try(read_header(fp), silent = TRUE)
  if (inherits(h, "try-error") || is.null(h)) next
  hm <- harmonise(h$chars)

  agents <- table(hm$Agent[!is.na(hm$Agent) & hm$Agent != "control"])
  nctrl  <- sum(!is.na(hm$Agent) & hm$Agent == "control")
  ok_ag  <- names(agents)[agents >= 2]

  rows[[length(rows) + 1]] <- data.frame(
    GSE = h$gse, Platform = h$gpl, Organism = h$org, nSamples = h$n,
    hasTreatmentField = any(!is.na(hm$TreatmentRaw)),
    nControl = nctrl, nAgents = length(ok_ag),
    Agents = paste(head(ok_ag, 6), collapse = "; "),
    hasCellLine = any(!is.na(hm$CellLine)),
    hasGenotype = any(!is.na(hm$Genotype)),
    usable = nctrl >= 2 && length(ok_ag) >= 1,
    Title = substr(h$title %||% "", 1, 90),
    stringsAsFactors = FALSE)
}
r <- do.call(rbind, rows)
dir.create("results", showWarnings = FALSE)
write.csv(r, file.path("results", "corpus_scan.csv"), row.names = FALSE)

cat(sprintf("\n=== Corpus scan: %d series in the bronze layer ===\n\n", nrow(r)))
cat(sprintf("  with any treatment field        %4d  (%.0f%%)\n",
            sum(r$hasTreatmentField), 100*mean(r$hasTreatmentField)))
cat(sprintf("  with a control arm (n>=2)       %4d  (%.0f%%)\n",
            sum(r$nControl >= 2), 100*mean(r$nControl >= 2)))
cat(sprintf("  USABLE: control + >=1 agent     %4d  (%.0f%%)\n",
            sum(r$usable), 100*mean(r$usable)))
cat(sprintf("  usable AND stratifiable by line %4d\n",
            sum(r$usable & (r$hasCellLine | r$hasGenotype))))
cat(sprintf("\n  total samples in usable studies %5d\n", sum(r$nSamples[r$usable])))
cat(sprintf("  distinct agents named           %5d\n",
            length(unique(unlist(strsplit(r$Agents[r$usable], "; "))))))

top <- sort(table(unlist(strsplit(r$Agents[r$usable], "; "))), decreasing = TRUE)
cat("\n  most frequently tested agents:\n")
for (i in seq_len(min(15, length(top))))
  cat(sprintf("    %-28s %d studies\n", names(top)[i], top[i]))

cat("\nWrote results/corpus_scan.csv\n")
