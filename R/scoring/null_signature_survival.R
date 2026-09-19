#!/usr/bin/env Rscript
# Is pLSC6 specifically prognostic, or would any six genes do?
#
# Venet, Dampier & Delattre (PLoS Comput Biol 2011) showed most random expression
# signatures are significantly associated with breast cancer outcome, because
# proliferation dominates the transcriptome. The same question has to be asked of any
# signature validated on survival. This compares pLSC6's prognostic performance in
# TARGET-AML against random six-gene signatures carrying the same coefficients.
#
#   Rscript R/null_signature_survival.R <bundle_dir> <manifest.csv> <clinical.csv> [nperm]

set.seed(20260917)
suppressPackageStartupMessages({ library(survival); library(data.table) })
source(file.path("R", "lsc_scores.R"))

a        <- commandArgs(trailingOnly = TRUE)
bundle   <- a[1]; manifest <- a[2]; clin_fp <- a[3]
NPERM    <- if (length(a) >= 4) as.integer(a[4]) else 1000L
CACHE    <- file.path("data", "processed", "target_gene_matrix.rds")

sig    <- load_signature("pLSC6")
genes  <- sig$gene
coefs  <- sig$coefficient
wanted <- unique(unlist(lapply(seq_len(nrow(sig)), function(i) symbols_for(sig[i, ]))))

if (file.exists(CACHE)) {
  M <- readRDS(CACHE)
  message(sprintf("Reusing cached matrix: %d genes x %d samples", nrow(M), ncol(M)))
} else {
  files <- list.files(bundle, pattern = "star_gene_counts\\.tsv$", recursive = TRUE,
                      full.names = TRUE)
  message(sprintf("Reading %d expression files ...", length(files)))
  lst <- lapply(files, function(f) {
    x <- fread(f, skip = 1, select = c("gene_name", "fpkm_unstranded"),
               showProgress = FALSE)
    x <- x[!is.na(fpkm_unstranded)]
    x[, .(v = max(fpkm_unstranded)), by = gene_name]   # collapse duplicate symbols
  })
  gs <- Reduce(intersect, lapply(lst, `[[`, "gene_name"))
  M  <- vapply(lst, function(x) log2(x[match(gs, gene_name), v] + 1), numeric(length(gs)))
  rownames(M) <- gs
  colnames(M) <- basename(dirname(files))

  M <- M[stats::complete.cases(M), , drop = FALSE]
  dir.create(dirname(CACHE), recursive = TRUE, showWarnings = FALSE)
  saveRDS(M, CACHE)
  message(sprintf("Built matrix: %d genes x %d samples", nrow(M), ncol(M)))
}

# alias mapping, as the GDC annotation uses ADGRG1 for GPR56
for (i in seq_len(nrow(sig))) {
  syn <- symbols_for(sig[i, ]); hit <- which(rownames(M) %in% syn)
  if (length(hit) && !(sig$gene[i] %in% rownames(M))) rownames(M)[hit[1]] <- sig$gene[i]
}
stopifnot(all(genes %in% rownames(M)))

man  <- read.csv(manifest, stringsAsFactors = FALSE)
clin <- read.csv(clin_fp, stringsAsFactors = FALSE)
clin$time  <- ifelse(clin$vital_status == "Dead", clin$days_to_death, clin$days_to_last_followup)
clin$event <- as.integer(clin$vital_status == "Dead")
clin <- clin[!is.na(clin$time) & clin$time >= 0, ]

idx <- data.frame(file_id = colnames(M), col = seq_len(ncol(M)), stringsAsFactors = FALSE)
d   <- merge(merge(idx, man, by = "file_id"), clin[, c("submitter_id","time","event")],
             by.x = "case", by.y = "submitter_id")
message(sprintf("Patients with expression and survival: %d (%d deaths)\n",
                nrow(d), sum(d$event)))

# Standardise every score so hazard ratios are per SD and therefore comparable.
perf <- function(idx_rows) {
  s  <- as.numeric(t(M[idx_rows, d$col, drop = FALSE]) %*% coefs)
  if (stats::sd(s) == 0) return(c(z = 0, C = 0.5))
  fit <- coxph(Surv(d$time, d$event) ~ scale(s))
  c(z = unname(summary(fit)$coefficients[, "z"]), C = unname(summary(fit)$concordance[1]))
}

obs <- perf(match(genes, rownames(M)))
cat(sprintf("pLSC6 observed:  Cox z = %.2f   C-index = %.3f\n\n", obs["z"], obs["C"]))

pool <- which(apply(M[, d$col, drop = FALSE], 1, stats::sd) > 0)
message(sprintf("Permuting %d random six-gene signatures from %d expressed genes ...",
                NPERM, length(pool)))
nul <- t(vapply(seq_len(NPERM), function(i) perf(sample(pool, length(genes))), numeric(2)))

pz <- (1 + sum(nul[, "z"]   >= obs["z"]))  / (NPERM + 1)   # higher z = worse survival
pC <- (1 + sum(nul[, "C"]   >= obs["C"]))  / (NPERM + 1)
cat(sprintf("Random signatures, Cox z: median %.2f, 95th pct %.2f\n",
            median(nul[, "z"]), quantile(nul[, "z"], .95)))
cat(sprintf("Random signatures, C-index: median %.3f, 95th pct %.3f\n\n",
            median(nul[, "C"]), quantile(nul[, "C"], .95)))
cat(sprintf("  pLSC6 beats random on Cox z:    p = %.4f\n", pz))
cat(sprintf("  pLSC6 beats random on C-index:  p = %.4f\n", pC))
cat(sprintf("  random signatures reaching P<0.05 in Cox: %.1f%%\n",
            100 * mean(abs(nul[, "z"]) > 1.96)))

saveRDS(list(obs = obs, null = nul), file.path("results", "null_survival_test.rds"))
write.csv(data.frame(nul), file.path("results", "null_survival_test.csv"), row.names = FALSE)
cat("\nWrote results/null_survival_test.csv\n")
