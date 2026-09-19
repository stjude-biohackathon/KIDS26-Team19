#!/usr/bin/env Rscript
# Figure: rank-based enrichment recovers signal a weighted sum cannot, and the agents it
# identifies group by mechanism.
#
# (A) The method comparison. Data, contrasts and permutation design are identical across all
#     four bars; only the summarisation and the gene set change.
# (B) Every contrast that lowered the LSC18 score and survived correction, grouped by agent.

r <- readRDS(file.path("results", "singscore_tests.rds"))
d <- read.csv(file.path("results", "lsc18_singscore_ranking.csv"), stringsAsFactors = FALSE)

ascii <- function(x) trimws(iconv(x, "UTF-8", "ASCII", sub = ""))
clip  <- function(x, n) ifelse(nchar(x) > n, paste0(substr(x, 1, n - 2), ".."), x)

s <- d[d$q < .10 & d$est < 0, ]
s$Agent <- clip(ascii(s$Agent), 24); s$Line <- clip(ascii(s$Line), 17)
ag <- names(sort(tapply(s$est, s$Agent, min)))
s$af <- factor(s$Agent, levels = ag); s <- s[order(s$af, s$est), ]

## panel A quantities, computed not asserted
sets  <- c("pLSC6", "LSC17", "LSC18")
surv  <- c(0, vapply(sets, function(n) sum(p.adjust(r[[n]]$null_p, "BH") < .10), numeric(1)))
beat  <- c(18, vapply(sets, function(n) sum(r[[n]]$null_p < .05), numeric(1)))
expd  <- c(12.9, vapply(sets, function(n) .05 * length(r[[n]]$null_p), numeric(1)))
fold  <- beat / expd
labs  <- c("weighted sum", "rank-based", "rank-based", "rank-based")
sublb <- c("pLSC6, 6 genes", "pLSC6, 6", "LSC17, 17", "LSC18, 18")

INK <- "#16191F"; MUT <- "#6B7480"; GRID <- "#E4E8ED"
ARR <- "#14532D"; CNT <- "#1E5A8A"; WS <- "#A8202B"

draw <- function() {
  layout(matrix(1:2, 1, 2), widths = c(1, 1.5))
  par(family = "sans")

  ## ---- A ---------------------------------------------------------------
  par(mar = c(6.6, 5.0, 5.0, 1.6), xpd = NA)
  ymax <- max(surv) * 1.30
  b <- barplot(surv, col = c(WS, ARR, ARR, ARR), border = NA, ylim = c(0, ymax),
               axes = FALSE, names.arg = rep("", 4), space = .45)
  axis(2, col = INK, col.axis = INK, las = 1, lwd = .9, cex.axis = .95,
       at = pretty(c(0, max(surv))))
  mtext("contrasts surviving FDR correction", side = 2, line = 3.1, cex = .9, col = INK)
  segments(b - .45, 0, b + .45, 0, col = INK, lwd = .9)

  text(b, surv + ymax * .045, surv, cex = 1.15, font = 2,
       col = c(WS, ARR, ARR, ARR))
  text(b, surv + ymax * .115, sprintf("%.1f×", fold), cex = .8, col = MUT)
  text(b, -ymax * .075, labs,  cex = .88, col = INK,  adj = .5)
  text(b, -ymax * .145, sublb, cex = .8,  col = MUT, adj = .5)

  mtext("A", side = 3, line = 3.2, adj = -0.22, cex = 1.25, font = 2, col = INK)
  mtext("Same data, same contrasts, same null", side = 3, line = 3.2, adj = 0,
        cex = 1.0, font = 2, col = INK)
  mtext("254 contrasts, 75 studies. A weighted sum has no internal reference, so any\nbroad perturbation moves it. Ranking each sample against its own transcriptome\nremoves that. Global P < 0.0005 for all three rank-based sets; the multiplier above\neach bar is enrichment over the number expected by chance.",
        side = 3, line = 0.35, adj = 0, cex = .74, col = MUT)

  ## ---- B ---------------------------------------------------------------
  par(mar = c(5.2, 16.5, 5.0, 3.0), xpd = NA)
  y  <- seq_len(nrow(s)) + cumsum(c(0, diff(as.integer(s$af)) != 0)) * 0.95
  xr <- c(min(s$est) * 1.08, -min(s$est) * 0.17)
  plot(NA, xlim = xr, ylim = c(max(y) + 1.6, min(y) - 1.6), axes = FALSE, xlab = "", ylab = "")
  at <- pretty(c(min(s$est), 0), 5)
  segments(at, min(y) - .8, at, max(y) + .5, col = GRID, lwd = .6)
  rect(0, y - .33, s$est, y + .33, col = ifelse(s$src == "array", ARR, CNT), border = NA)
  segments(0, min(y) - .8, 0, max(y) + .5, col = INK, lwd = 1.1)

  axis(1, at = at, col = INK, col.axis = INK, lwd = .9, cex.axis = .92,
       pos = max(y) + 1.0)
  mtext("shift in LSC18 enrichment vs control (singscore units)", side = 1, line = 2.2,
        cex = .9, col = INK)

  text(xr[1] * 1.03, y, s$Line, adj = 1, cex = .84, col = INK)
  text(xr[2] * 0.14, y, sprintf("q %.3f", s$q), adj = 0, cex = .78, col = MUT)
  for (a in levels(s$af)) {
    i <- which(s$af == a)
    text(xr[1] * 1.52, mean(y[i]), a, adj = 1, cex = .9, font = 2, col = INK)
    if (length(i) > 1)
      lines(rep(xr[1] * 1.42, 2), c(min(y[i]) - .42, max(y[i]) + .42),
            col = "#C3CBD8", lwd = 2.2)
  }

  mtext("B", side = 3, line = 3.2, adj = -0.30, cex = 1.25, font = 2, col = INK)
  mtext("Agents lowering the stemness score", side = 3, line = 3.2, adj = 0,
        cex = 1.0, font = 2, col = INK)
  nb <- sum(s$Agent == "BAY-1251152")
  mtext(sprintf("All %d contrasts surviving FDR in this direction. BAY-1251152, a menin-MLL\ninhibitor, lowers the score in %d cell lines; four carry KMT2A rearrangements.\nEleven further contrasts survived in the opposite direction and are not shown.",
                nrow(s), nb),
        side = 3, line = 0.35, adj = 0, cex = .74, col = MUT)

  ## legend above the plotting band. Offsets are absolute, because xr[1] is negative and
  ## scaling by it flips the direction the labels move.
  ly  <- min(y) - 1.05
  w   <- abs(xr[1]) * 0.035
  gap <- abs(xr[1]) * 0.022
  lx  <- xr[1] * 1.52
  rect(lx, ly - .26, lx + w, ly + .26, col = ARR, border = NA)
  text(lx + w + gap, ly, "array series matrix", adj = 0, cex = .82, col = INK)
  lx2 <- xr[1] * 0.82
  rect(lx2, ly - .26, lx2 + w, ly + .26, col = CNT, border = NA)
  text(lx2 + w + gap, ly, "NCBI RNA-seq counts", adj = 0, cex = .82, col = INK)
}

H <- max(6.8, 3.9 + nrow(s) * 0.30)
pdf(file.path("results", "figure_singscore.pdf"), width = 13.6, height = H, useDingbats = FALSE)
draw(); invisible(dev.off())
png(file.path("results", "figure_singscore.png"), width = 13.6, height = H, units = "in", res = 600)
draw(); invisible(dev.off())
cat(sprintf("panel A: %s survivors | panel B: %d contrasts, %d agents\n",
            paste(surv, collapse = "/"), nrow(s), nlevels(s$af)))
