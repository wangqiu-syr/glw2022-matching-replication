# plots.R -- Stage 3 figures and the Simulation 2 comparison with the paper.
# Base R graphics only. Run from the repository root after R/summarize.R:
#   Rscript R/plots.R
#
# Inputs (results/):
#   sim2_summary.csv      from R/summarize.R
#   paper_table3.csv      paper's Table 3 (Simulation 2), published values. Columns:
#                           pair ("t,tp"), DIF_bias, DIF_mse, CM_bias, CM_mse, CM_coverage,
#                           GPS_bias, GPS_mse, GPS_coverage   [DIF_coverage optional; the paper has none]
#                         (bias = signed bias; coverage in %). Row order = Fig. 1 index.
#   paper_table5.csv      paper's Table 5 (application, 15 effects), published values. Columns:
#                           pair, DIF_est, CM_est, CM_LL, CM_UL, GPS_est, GPS_LL, GPS_UL
#                           [CM_sd, GPS_sd present but not used]
# Outputs (results/):
#   fig1_bias_sim2_pop.png, fig1_bias_sim2_bal.png      Fig. 1 replica (ours)
#   fig1_bias_sim2_pop_vs_paper.png                     ours + paper overlay  (needs table 3)
#   fig2_table5_redrawn.png                             Fig. 2 redrawn        (needs table 5)
#   sim2_vs_paper.csv                                   per pair x estimator  (needs table 3)
# Parts whose paper table is missing are skipped with a message.
# Optional argument: a different folder for inputs and outputs (default results).

RES <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) a[1] else "results" }
EST <- c("DIF", "CM", "GPS")
COL <- c(DIF = "red", CM = "green3", GPS = "blue")
LWD_PAPER <- 1.3   # width of the paper's dashed lines in the overlays (plot and legend)

read_paper <- function(file, required) {
  f <- file.path(RES, file)
  if (!file.exists(f)) {
    cat(sprintf("SKIP: %s not found. Expected columns: %s\n", f, paste(required, collapse = ", ")))
    return(NULL)
  }
  x <- read.csv(f, check.names = FALSE)
  miss <- setdiff(required, names(x))
  if (length(miss)) {
    cat(sprintf("SKIP: %s lacks columns: %s\n", f, paste(miss, collapse = ", ")))
    return(NULL)
  }
  x
}

summ <- read.csv(file.path(RES, "sim2_summary.csv"))
summ$pair <- paste0(summ$t, ",", summ$tp)

t3_cols <- c("pair", as.vector(outer(EST, c("_bias", "_mse"), paste0)),
             "CM_coverage", "GPS_coverage")
t3 <- read_paper("paper_table3.csv", t3_cols)

# Index order for Fig. 1: the paper table's row order if available, else (1,2), (1,3), ..., (8,9).
pairs <- if (!is.null(t3)) gsub(" ", "", t3$pair) else {
  pr <- t(combn(9, 2)); paste0(pr[, 1], ",", pr[, 2])
}
stopifnot(length(pairs) == 36, setequal(pairs, unique(summ$pair)))

ours_bias <- function(target) {
  sapply(EST, function(e) {
    x <- summ[summ$estimator == e & summ$target == target, ]
    x$bias[match(pairs, x$pair)]
  })
}

fig1 <- function(file, target, paper = NULL, paper_label = "paper", subtitle = NULL) {
  b <- ours_bias(target)
  yl <- range(c(0, b, if (!is.null(paper)) paper), na.rm = TRUE)
  png(file.path(RES, file), width = 1600, height = 900, res = 150)
  op <- par(mar = c(5.5, 4.5, 3.5, 1))
  plot(NA, xlim = c(1, 36), ylim = yl, xaxt = "n", xlab = "",
       ylab = sprintf("Signed bias (vs tau_%s)", target),
       main = sprintf("Simulation 2: bias of the 36 pairwise effects (500 reps, target tau_%s)", target))
  axis(1, at = 1:36, labels = paste0("(", pairs, ")"), las = 2, cex.axis = 0.65)
  mtext("Effect index (pair t,t')", side = 1, line = 4.3)
  if (!is.null(subtitle)) mtext(subtitle, side = 3, line = 0.4, cex = 0.85)
  abline(h = 0, col = "black", lwd = 1.5)
  for (e in EST) lines(1:36, b[, e], type = "o", pch = 16, cex = 0.7, lwd = 1.8, col = COL[[e]])
  if (!is.null(paper)) {
    for (e in EST) lines(1:36, paper[, e], lty = 2, lwd = LWD_PAPER, col = COL[[e]])
    legend("topleft", bty = "n", ncol = 2, cex = 0.85,
           legend = c(paste(EST, "ours"), paste(EST, paper_label)),
           col = unname(c(COL[EST], COL[EST])), lty = rep(c(1, 2), each = 3),
           lwd = rep(c(1.8, LWD_PAPER), each = 3), seg.len = 2.5,
           pch = rep(c(16, NA), each = 3))
  } else {
    legend("topleft", bty = "n", legend = EST, col = COL[EST], lty = 1, lwd = 1.8, pch = 16, cex = 0.85)
  }
  par(op); dev.off()
  cat(sprintf("wrote %s\n", file.path(RES, file)))
}

# A1. Fig. 1 replica, both targets.
fig1("fig1_bias_sim2_pop.png", "pop")
fig1("fig1_bias_sim2_bal.png", "bal")

if (!is.null(t3)) {
  # A2. Overlay with the paper's values. The paper's truth is the population
  # ATE (all 0 in Simulation 2), so the overlay uses tau_pop.
  paper_b <- sapply(EST, function(e) t3[[paste0(e, "_bias")]])
  fig1("fig1_bias_sim2_pop_vs_paper.png", "pop", paper = paper_b,
       subtitle = "paper's bias as published")
  # Same overlay with the paper's bias multiplied by -1 (see docs/spec.md,
  # "Comparison with the paper"). The published-sign version above is kept.
  fig1("fig1_bias_sim2_pop_vs_paper_flipped.png", "pop", paper = -paper_b,
       paper_label = "paper x (-1)", subtitle = "paper's bias with its sign flipped (x -1)")

  # Correlation across the 36 pairs between our bias and the paper's bias,
  # with the paper's sign as published and flipped (flipping only negates r).
  cors <- do.call(rbind, lapply(c("pop", "bal"), function(target) {
    b <- ours_bias(target)
    do.call(rbind, lapply(EST, function(e) data.frame(
      target = target, estimator = e,
      r_published = cor(b[, e], paper_b[, e]),
      r_flipped   = cor(b[, e], -paper_b[, e]),
      same_sign_published = sum(sign(b[, e]) == sign(paper_b[, e])),
      same_sign_flipped   = sum(sign(b[, e]) == sign(-paper_b[, e])),
      n_pairs = 36)))
  }))
  write.csv(cors, file.path(RES, "sim2_bias_correlations.csv"), row.names = FALSE)
  cat("\nCorrelation of our bias with the paper's bias across the 36 pairs:\n")
  print(format(cors, digits = 3), row.names = FALSE)
  cat(sprintf("wrote %s\n", file.path(RES, "sim2_bias_correlations.csv")))

  # A4. Pair-by-pair comparison.
  cmp <- do.call(rbind, lapply(EST, function(e) {
    get <- function(target, col) {
      x <- summ[summ$estimator == e & summ$target == target, ]
      x[[col]][match(pairs, x$pair)]
    }
    data.frame(index = 1:36, pair = paste0("(", pairs, ")"), estimator = e,
               paper_bias = t3[[paste0(e, "_bias")]],
               ours_bias_pop = get("pop", "bias"), ours_bias_bal = get("bal", "bias"),
               paper_mse = t3[[paste0(e, "_mse")]],
               ours_mse_pop = get("pop", "mse"), ours_mse_bal = get("bal", "mse"),
               paper_coverage = if (paste0(e, "_coverage") %in% names(t3))
                                  t3[[paste0(e, "_coverage")]] else NA_real_,
               ours_coverage_pop = get("pop", "coverage"), ours_coverage_bal = get("bal", "coverage"),
               n_reps = get("pop", "n_reps"))
  }))
  write.csv(cmp, file.path(RES, "sim2_vs_paper.csv"), row.names = FALSE)
  cat(sprintf("wrote %s\n", file.path(RES, "sim2_vs_paper.csv")))
}

# A3. Fig. 2 redrawn from the paper's Table 5 (not a replication).
t5 <- read_paper("paper_table5.csv",
                 c("pair", "DIF_est", "CM_est", "CM_LL", "CM_UL", "GPS_est", "GPS_LL", "GPS_UL"))
if (!is.null(t5)) {
  k <- nrow(t5)
  png(file.path(RES, "fig2_table5_redrawn.png"), width = 1600, height = 900, res = 150)
  op <- par(mar = c(6, 4.5, 3.5, 1))
  yl <- range(t5[, c("DIF_est", "CM_est", "CM_LL", "CM_UL", "GPS_est", "GPS_LL", "GPS_UL")], na.rm = TRUE)
  plot(NA, xlim = c(1, k), ylim = yl, xaxt = "n", xlab = "", ylab = "Estimated effect",
       main = "Redrawn from the paper's Table 5 (not replicated)")
  axis(1, at = 1:k, labels = paste0("(", gsub(" ", "", t5$pair), ")"), las = 2, cex.axis = 0.75)
  mtext("Treatment pair", side = 1, line = 4.3)
  abline(h = 0, col = "grey50")
  lines(1:k, t5$DIF_est, type = "o", pch = 16, lwd = 1.8, col = COL[["DIF"]])
  lines(1:k, t5$CM_est,  type = "o", pch = 16, lwd = 1.8, col = COL[["CM"]])
  lines(1:k, t5$GPS_est, type = "o", pch = 16, lwd = 1.8, col = COL[["GPS"]])
  for (b in c("CM_LL", "CM_UL"))   lines(1:k, t5[[b]], lty = 2, col = COL[["CM"]])
  for (b in c("GPS_LL", "GPS_UL")) lines(1:k, t5[[b]], lty = 2, col = COL[["GPS"]])
  legend("topleft", bty = "n", cex = 0.85,
         legend = c("DIF estimate", "CM estimate", "GPS estimate", "CM 95% limits", "GPS 95% limits"),
         col = COL[c("DIF", "CM", "GPS", "CM", "GPS")], lty = c(1, 1, 1, 2, 2), pch = c(16, 16, 16, NA, NA))
  par(op); dev.off()
  cat(sprintf("wrote %s\n", file.path(RES, "fig2_table5_redrawn.png")))
}
