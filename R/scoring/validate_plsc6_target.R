#!/usr/bin/env Rscript
# Does our pLSC6 implementation reproduce the published TARGET validation?
#
# Elsayed et al. (Leukemia 2020, doi:10.1038/s41375-019-0604-8) validated pLSC6 in
# 205 TARGET pediatric AML patients and reported, for the high vs low score groups:
#   OS  HR = 2.81 (95% CI 1.85-4.28)
#   EFS HR = 2.86 (95% CI 2.02-4.04); 5-yr EFS 49.2% low vs 13.7% high
#   split by recursive partitioning into 60% low / 40% high
#
# Only overall survival is reproduced here: EFS lives in the TARGET clinical
# supplements, not in the case-level GDC fields. Matching their direction and rough
# magnitude is evidence the scoring code is correct; it is not a new finding.
#
#   Rscript R/validate_plsc6_target.R <bundle_dir> <manifest.csv> <clinical.csv>

suppressPackageStartupMessages({ library(survival) })
source(file.path("R", "lsc_scores.R"))

args     <- commandArgs(trailingOnly = TRUE)
bundle   <- if (length(args) >= 1) args[1] else "target/bundle"
manifest <- if (length(args) >= 2) args[2] else "target/manifest.csv"
clin_fp  <- if (length(args) >= 3) args[3] else "data/processed/target_aml_clinical.csv"

sig   <- load_signature("pLSC6")
stopifnot(all(sig$verified))
wanted <- unique(unlist(lapply(seq_len(nrow(sig)), function(i) symbols_for(sig[i, ]))))

man <- read.csv(manifest, stringsAsFactors = FALSE)
cat(sprintf("manifest: %d cases\n", nrow(man)))

# Each GDC STAR file is one sample: gene_id, gene_name, ..., tpm, fpkm.
# Elsayed used log2(RPKM + 1); fpkm_unstranded is the closest GDC equivalent.
read_star <- function(fp) {
  x <- read.delim(fp, comment.char = "#", stringsAsFactors = FALSE, check.names = FALSE)
  x <- x[!grepl("^N_", x$gene_id), , drop = FALSE]
  hit <- x[x$gene_name %in% wanted, c("gene_name", "fpkm_unstranded"), drop = FALSE]
  if (!nrow(hit)) return(NULL)
  v <- tapply(hit$fpkm_unstranded, hit$gene_name, max, na.rm = TRUE)
  log2(v + 1)
}

files <- list.files(bundle, pattern = "star_gene_counts\\.tsv$", recursive = TRUE, full.names = TRUE)
cat(sprintf("expression files found: %d\n", length(files)))
if (!length(files)) stop("No STAR count files under ", bundle, call. = FALSE)

uuid_of <- function(p) basename(dirname(p))
rows <- list()
for (f in files) {
  v <- read_star(f)
  if (is.null(v)) next
  rows[[uuid_of(f)]] <- v
}
genes <- Reduce(intersect, lapply(rows, names))
cat(sprintf("samples parsed: %d;  signature genes present in all: %d/%d (%s)\n",
            length(rows), length(genes), nrow(sig), paste(genes, collapse = ", ")))

# Map whichever alias the annotation used back to the signature's primary name.
owner <- vapply(genes, function(s) {
  hit <- which(vapply(seq_len(nrow(sig)), function(i) s %in% symbols_for(sig[i, ]), logical(1)))
  if (length(hit)) sig$gene[hit[1]] else NA_character_
}, character(1))
genes <- genes[!is.na(owner)]; owner <- owner[!is.na(owner)]
if (length(genes) < nrow(sig)) {
  stop(sprintf("Only %d of %d pLSC6 genes available; refusing to score a partial signature.",
               length(genes), nrow(sig)), call. = FALSE)
}

mat <- t(vapply(rows, function(v) v[genes], numeric(length(genes))))
colnames(mat) <- owner
coefs <- setNames(sig$coefficient, sig$gene)[colnames(mat)]
score <- as.numeric(mat %*% coefs)

df <- data.frame(file_id = rownames(mat), pLSC6 = score, stringsAsFactors = FALSE)
df <- merge(df, man, by = "file_id")

clin <- read.csv(clin_fp, stringsAsFactors = FALSE)
clin$time  <- ifelse(clin$vital_status == "Dead", clin$days_to_death, clin$days_to_last_followup)
clin$event <- as.integer(clin$vital_status == "Dead")
clin <- clin[!is.na(clin$time) & !is.na(clin$event) & clin$time >= 0, ]

d <- merge(df, clin[, c("submitter_id", "time", "event")],
           by.x = "case", by.y = "submitter_id")
cat(sprintf("scored cases joined to survival: %d (%d deaths)\n\n", nrow(d), sum(d$event)))

# The paper dichotomised by recursive partitioning and landed on 60% low / 40% high.
# Applying that split directly keeps the comparison like-for-like.
cut60 <- quantile(d$pLSC6, 0.60, na.rm = TRUE)
d$group <- factor(ifelse(d$pLSC6 > cut60, "high", "low"), levels = c("low", "high"))

cat(sprintf("pLSC6: median %.2f, range %.2f..%.2f; cut at 60th pct = %.2f\n",
            median(d$pLSC6), min(d$pLSC6), max(d$pLSC6), cut60))
print(table(d$group))

fit <- survfit(Surv(time / 365.25, event) ~ group, data = d)
s5  <- summary(fit, times = 5, extend = TRUE)
cat("\n5-year overall survival\n")
for (i in seq_along(s5$strata)) {
  cat(sprintf("  %-12s %.3f (95%% CI %.3f-%.3f)\n",
              sub("group=", "", as.character(s5$strata[i])), s5$surv[i], s5$lower[i], s5$upper[i]))
}

cox <- coxph(Surv(time, event) ~ group, data = d)
hr  <- summary(cox)
cat(sprintf("\nOS hazard ratio, high vs low: %.2f (95%% CI %.2f-%.2f), p = %.3g\n",
            hr$conf.int[1], hr$conf.int[3], hr$conf.int[4], hr$coefficients[5]))

cont <- coxph(Surv(time, event) ~ pLSC6, data = d)
hc   <- summary(cont)
cat(sprintf("per-unit pLSC6 (continuous):   %.2f (95%% CI %.2f-%.2f), p = %.3g\n",
            hc$conf.int[1], hc$conf.int[3], hc$conf.int[4], hc$coefficients[5]))

cat(sprintf("concordance (C-index):         %.3f\n", hc$concordance[1]))
cat("\nPublished TARGET validation (Elsayed 2020): OS HR = 2.81 (1.85-4.28), n = 205\n")

dir.create("results", showWarnings = FALSE)

# Vector for print, raster for the README and slides.
km_plot <- function() {
  op <- par(mar = c(4.6, 4.8, 3.4, 1.6), family = "sans")
  plot(fit, col = c("#2B6CB0", "#9B2C2C"), lwd = 2.4, mark.time = TRUE,
       xlab = "", ylab = "", axes = FALSE)
  axis(1, col = "#16191F", col.axis = "#16191F", lwd = .9, cex.axis = .95)
  axis(2, col = "#16191F", col.axis = "#16191F", lwd = .9, las = 1, cex.axis = .95)
  mtext("Years from diagnosis", side = 1, line = 2.9, cex = 1.0, col = "#16191F")
  mtext("Overall survival", side = 2, line = 3.2, cex = 1.0, col = "#16191F")
  mtext(sprintf("TARGET-AML overall survival by pLSC6 group (n = %d)", nrow(d)),
        side = 3, line = 1.7, adj = 0, cex = 1.15, font = 2, col = "#16191F")
  mtext(sprintf("%d deaths. Split at the 60th percentile, matching the source publication.",
                sum(d$event)),
        side = 3, line = 0.4, adj = 0, cex = .82, col = "#6B7480")
  legend("bottomleft", bty = "n", cex = .92, lwd = 2.4, col = c("#2B6CB0", "#9B2C2C"),
         legend = c(sprintf("low pLSC6, 60%% (n = %d)",  sum(d$group == "low")),
                    sprintf("high pLSC6, 40%% (n = %d)", sum(d$group == "high"))))
  par(op)
}

pdf(file.path("results", "plsc6_target_km.pdf"), width = 7.2, height = 5.2,
    useDingbats = FALSE); km_plot(); invisible(dev.off())
png(file.path("results", "plsc6_target_km.png"), width = 7.2, height = 5.2,
    units = "in", res = 600); km_plot(); invisible(dev.off())
cat("\nWrote results/plsc6_target_km.{pdf,png}\n")
