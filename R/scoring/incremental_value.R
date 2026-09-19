#!/usr/bin/env Rscript
# Does pLSC6 add prognostic information beyond what a clinician already has?
#
# A score that merely restates age or sex adds nothing in practice. This fits nested
# Cox models on the TARGET-AML cohort and asks whether pLSC6 survives adjustment, and
# whether adding it to a clinical baseline improves discrimination.
#
#   Rscript R/incremental_value.R <manifest.csv> <clinical.csv>

suppressPackageStartupMessages(library(survival))
source(file.path("R", "lsc_scores.R"))

a    <- commandArgs(trailingOnly = TRUE)
man  <- read.csv(a[1], stringsAsFactors = FALSE)
clin <- read.csv(a[2], stringsAsFactors = FALSE)

sig <- load_signature("pLSC6")
M   <- readRDS(file.path("data", "processed", "target_gene_matrix.rds"))
for (i in seq_len(nrow(sig))) {
  syn <- symbols_for(sig[i, ]); h <- which(rownames(M) %in% syn)
  if (length(h) && !(sig$gene[i] %in% rownames(M))) rownames(M)[h[1]] <- sig$gene[i]
}
score <- as.numeric(t(M[sig$gene, , drop = FALSE]) %*% sig$coefficient)

clin$time  <- ifelse(clin$vital_status == "Dead", clin$days_to_death, clin$days_to_last_followup)
clin$event <- as.integer(clin$vital_status == "Dead")
clin <- clin[!is.na(clin$time) & clin$time >= 0, ]

d <- merge(data.frame(file_id = colnames(M), pLSC6 = score, stringsAsFactors = FALSE),
           man, by = "file_id")
d <- merge(d, clin[, c("submitter_id","time","event","age_at_diagnosis_days","gender")],
           by.x = "case", by.y = "submitter_id")
d$age   <- d$age_at_diagnosis_days / 365.25
d$pLSC6 <- as.numeric(scale(d$pLSC6))
d <- d[!is.na(d$age), ]
cat(sprintf("Cohort: %d patients, %d deaths, age %.1f-%.1f years (median %.1f)\n\n",
            nrow(d), sum(d$event), min(d$age), max(d$age), median(d$age)))

cat(sprintf("Correlation of pLSC6 with age: r = %.3f\n\n",
            cor(d$pLSC6, d$age, use = "complete.obs")))

S <- Surv(d$time, d$event)
m_age  <- coxph(S ~ age, data = d)
m_lsc  <- coxph(S ~ pLSC6, data = d)
m_both <- coxph(S ~ age + pLSC6, data = d)

row <- function(fit, term, label) {
  s <- summary(fit); i <- match(term, rownames(s$coefficients))
  cat(sprintf("  %-34s HR %.2f (%.2f-%.2f)  P = %-8.2g\n", label,
              s$conf.int[i,1], s$conf.int[i,3], s$conf.int[i,4], s$coefficients[i,5]))
}
cat("Single-predictor models\n")
row(m_age, "age",   "age (per year)")
row(m_lsc, "pLSC6", "pLSC6 (per SD)")
cat("\nAdjusted model: age + pLSC6\n")
row(m_both, "age",   "age (per year)")
row(m_both, "pLSC6", "pLSC6 (per SD), adjusted for age")

lr <- anova(m_age, m_both)
cat(sprintf("\nLikelihood ratio, adding pLSC6 to age: chisq = %.1f, df = %d, P = %.3g\n",
            lr$Chisq[2], lr$Df[2], lr$`Pr(>|Chi|)`[2]))
lr2 <- anova(m_lsc, m_both)
cat(sprintf("Likelihood ratio, adding age to pLSC6: chisq = %.1f, df = %d, P = %.3g\n",
            lr2$Chisq[2], lr2$Df[2], lr2$`Pr(>|Chi|)`[2]))

cat("\nDiscrimination (C-index)\n")
for (nm in c("age","pLSC6","age + pLSC6")) {
  f <- switch(nm, "age" = m_age, "pLSC6" = m_lsc, m_both)
  cat(sprintf("  %-14s %.3f\n", nm, summary(f)$concordance[1]))
}
