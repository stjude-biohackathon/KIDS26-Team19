#!/usr/bin/env Rscript
# Where the signature carries specific information, and where it does not.
#
# Three permutation tests against random six-gene signatures carrying pLSC6's published
# coefficients. Panels B and C are independent corpora: array-based series matrices and
# NCBI-recomputed RNA-seq counts, differing in platform, quantification and normalisation.
# They agree, which is what makes the negative result a replication rather than a repetition.

sv <- readRDS(file.path("results", "null_survival_test.rds"))
ar <- readRDS(file.path("results", "corpus_null_draws.rds"))
nc <- readRDS(file.path("results", "ncbi_null_draws.rds"))

INK<-"#16191F"; MUT<-"#6B7480"; GRID<-"#E4E8ED"
NULLC<-"#C7D0DA"; NULLE<-"#9AA6B2"; HIT<-"#14532D"; MISS<-"#A8202B"

panel <- function(null, obs, side, xlab, title, sub, lab, pval, hit) {
  h <- hist(null, breaks = 34, plot = FALSE)
  xr <- range(c(h$breaks, obs)) + c(-.05, .05) * diff(range(c(h$breaks, obs)))
  plot(NA, xlim = xr, ylim = c(0, max(h$counts) * 1.32), axes = FALSE, xlab = "", ylab = "")
  sel <- if (side == "right") h$mids >= obs else h$mids <= obs
  rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$counts,
       col = ifelse(sel, if (hit) HIT else MISS, NULLC), border = NULLE, lwd = .4)
  axis(1, col = INK, col.axis = INK, cex.axis = .95, lwd = .8, tck = -.022, mgp = c(3,.7,0))
  axis(2, col = INK, col.axis = INK, cex.axis = .95, lwd = .8, las = 1, tck = -.022, mgp = c(3,.6,0))
  mtext(xlab, side = 1, line = 2.5, cex = .82, col = INK)
  mtext("random signatures", side = 2, line = 2.9, cex = .82, col = INK)
  segments(obs, 0, obs, max(h$counts)*1.11, col = if (hit) HIT else MISS, lwd = 2.4)
  points(obs, max(h$counts)*1.11, pch = 25, bg = if (hit) HIT else MISS,
         col = if (hit) HIT else MISS, cex = 1.15)
  adj <- if (obs > mean(xr)) 1.03 else -0.03
  text(obs, max(h$counts)*1.21, lab, col = if (hit) HIT else MISS, cex = .95, font = 2, adj = adj)
  text(obs, max(h$counts)*1.21, sprintf("  P = %s  ", pval), col = INK, cex = .88,
       adj = c(adj, 1.9))
  mtext(title, side = 3, line = 1.25, adj = 0, cex = .98, font = 2, col = INK)
  mtext(sub,   side = 3, line = 0.1,  adj = 0, cex = .72, col = MUT)
}

pv <- function(null, obs, side = "right") {
  v <- (1 + sum(if (side == "right") null >= obs else null <= obs)) / (length(null) + 1)
  if (v < 0.01) sprintf("%.3f", v) else sprintf("%.2f", v)
}

draw <- function() {
  par(mfrow = c(1,3), mar = c(4.4,4.6,4.2,1.4), family = "sans", oma = c(0.4,0.4,2.6,0.4))
  panel(sv$null[,"z"], sv$obs["z"], "right", "Cox z statistic",
        "A   Survival, 351 TARGET patients",
        "31% of random signatures reach P < 0.05 here",
        "pLSC6 4.43", pv(sv$null[,"z"], sv$obs["z"]), TRUE)
  panel(ar$global_null, ar$global_obs, "right", "mean absolute effect, 95 contrasts",
        "B   Drug response, array series",
        "38 series, 41 agents: random sets do better",
        sprintf("pLSC6 %.2f", ar$global_obs), pv(ar$global_null, ar$global_obs), FALSE)
  panel(nc$global_null, nc$global_obs, "right", "mean absolute effect, 162 contrasts",
        "C   Drug response, NCBI RNA-seq counts",
        "38 further series, 74 agents: the same answer",
        sprintf("pLSC6 %.2f", nc$global_obs), pv(nc$global_null, nc$global_obs), FALSE)
  mtext("pLSC6 against 1,000 random six-gene signatures carrying the same coefficients. Panels B and C are independent corpora.",
        side = 3, outer = TRUE, line = 0.5, adj = 0, cex = .92, col = INK, font = 2)
}

pdf(file.path("results","figure_null_tests.pdf"), width = 12.6, height = 4.5, useDingbats = FALSE)
draw(); invisible(dev.off())
png(file.path("results","figure_null_tests.png"), width = 12.6, height = 4.5, units = "in", res = 600)
draw(); invisible(dev.off())
cat(sprintf("A: %.2f (P=%s)  B: %.3f vs %.3f (P=%s)  C: %.3f vs %.3f (P=%s)\n",
  sv$obs["z"], pv(sv$null[,"z"], sv$obs["z"]),
  ar$global_obs, median(ar$global_null), pv(ar$global_null, ar$global_obs),
  nc$global_obs, median(nc$global_null), pv(nc$global_null, nc$global_obs)))
