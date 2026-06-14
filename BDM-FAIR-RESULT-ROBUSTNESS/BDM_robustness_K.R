# BDM_robustness_K_translated.R
# Translation of BDM_robustness_K.m from MATLAB to R
# Generates Figures 7 and 8 using the uploaded .mat result workspaces.

# -----------------------------
# Setup
# -----------------------------
rm(list = ls())

required_packages <- c("R.matlab")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(R.matlab)

# Change this if the .mat files are stored elsewhere.
DATA_DIR <- "."
OUTPUT_DIR <- file.path(DATA_DIR, "R_OUTPUT")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# -----------------------------
# Helper functions
# -----------------------------
load_msumS <- function(file_name) {
  file_path <- file.path(DATA_DIR, file_name)
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  mat <- R.matlab::readMat(file_path)

  if (!"MsumS" %in% names(mat)) {
    stop("MsumS not found in: ", file_name)
  }

  mat$MsumS
}

compute_results_prob <- function(MsumS, nhor) {
  # MATLAB equivalent:
  # count = sum(squeeze(MsumS(:,2,hor_Msum)) > squeeze(MsumS(:,1,hor_Msum)));
  # results_prob(hor_Msum) = 1 - count / length(MsumS);

  n_draws <- dim(MsumS)[1]
  results_prob <- numeric(nhor)

  for (hor_Msum in seq_len(nhor)) {
    count <- sum(MsumS[, 2, hor_Msum] > MsumS[, 1, hor_Msum])
    results_prob[hor_Msum] <- 1 - count / n_draws
  }

  results_prob
}

# -----------------------------
# BDM baseline output
# -----------------------------
nhor <- 30

orig <- compute_results_prob(
  load_msumS("asym_results_BDM.mat"),
  nhor = 30
)

orig_30 <- compute_results_prob(
  load_msumS("asym_results_BDM_30lags.mat"),
  nhor = 30
)

orig_20 <- compute_results_prob(
  load_msumS("asym_results_BDM_20lags.mat"),
  nhor = 20
)

# -----------------------------
# Outputs from the full RZ sample, 1890-2015
# -----------------------------
origfullsample <- compute_results_prob(
  load_msumS("asym_results_fullRZsample.mat"),
  nhor = 30
)

orig_30fullsample <- compute_results_prob(
  load_msumS("asym_results_fullRZsample_30lags.mat"),
  nhor = 30
)

orig_20fullsample <- compute_results_prob(
  load_msumS("asym_results_fullRZsample_20lags.mat"),
  nhor = 20
)

# -----------------------------
# Figure 7
# -----------------------------
pdf(file.path(OUTPUT_DIR, "figure_7_BDM_FAIR_robustness.pdf"), width = 8, height = 6)

plot(
  seq_len(nhor), orig,
  type = "l",
  lwd = 1.5,
  xlab = "horizon, h",
  ylab = "",
  main = "P(msum neg>msum pos) at hor=h",
  ylim = range(c(orig, origfullsample, 0.9, 0.95), na.rm = TRUE)
)
lines(seq_len(nhor), origfullsample, lwd = 1.5)
lines(seq_len(nhor), rep(0.9, nhor), col = "black")
lines(seq_len(nhor), rep(0.95, nhor), col = "black")
legend(
  "bottomright",
  legend = c("BDM sample: 1901q2-2015q4", "RZ full sample: 1890q1-2015q4", "", ""),
  lty = c(1, 1, 1, 1),
  lwd = c(1.5, 1.5, 1, 1),
  col = c("black", "black", "black", "black"),
  bty = "n"
)

dev.off()

# -----------------------------
# Figure 8
# -----------------------------
pdf(file.path(OUTPUT_DIR, "figure_8_BDM_FAIR_robustness.pdf"), width = 12, height = 6)

par(mfrow = c(1, 2))

plot(
  seq_len(nhor), orig,
  type = "l",
  col = "blue",
  lwd = 2,
  xlab = "horizon, h",
  ylab = "",
  main = "1901q2-2015q4",
  ylim = range(c(orig, orig_30, orig_20, 0.9, 0.95), na.rm = TRUE)
)
lines(seq_len(nhor), orig_30, col = "red", lty = 2, lwd = 2)
lines(seq_along(orig_20), orig_20, col = "cyan", type = "b", pch = 8, lwd = 2)
lines(seq_len(nhor), rep(0.9, nhor), col = "black", lwd = 1)
lines(seq_len(nhor), rep(0.95, nhor), col = "black", lwd = 1)

plot(
  seq_len(nhor), origfullsample,
  type = "l",
  col = "blue",
  lwd = 2,
  xlab = "horizon, h",
  ylab = "",
  main = "1890q1-2015q4",
  ylim = range(c(origfullsample, orig_30fullsample, orig_20fullsample, 0.9, 0.95), na.rm = TRUE)
)
lines(seq_len(nhor), orig_30fullsample, col = "red", lty = 2, lwd = 2)
lines(seq_along(orig_20fullsample), orig_20fullsample, col = "cyan", type = "b", pch = 8, lwd = 2)
lines(seq_len(nhor), rep(0.9, nhor), col = "black", lwd = 1)
lines(seq_len(nhor), rep(0.95, nhor), col = "black", lwd = 1)
legend(
  "bottomright",
  legend = c("Original specification w/ 45lags", "w/ 30 lags", "w/ 20 lags", "", ""),
  col = c("blue", "red", "cyan", "black", "black"),
  lty = c(1, 2, 1, 1, 1),
  pch = c(NA, NA, 8, NA, NA),
  lwd = c(2, 2, 2, 1, 1),
  bty = "n"
)

dev.off()

par(mfrow = c(1, 1))

# Optional: export the computed probability series for checking.
results_df <- data.frame(
  horizon = seq_len(nhor),
  orig = orig,
  orig_30 = orig_30,
  orig_20 = c(orig_20, rep(NA_real_, nhor - length(orig_20))),
  origfullsample = origfullsample,
  orig_30fullsample = orig_30fullsample,
  orig_20fullsample = c(orig_20fullsample, rep(NA_real_, nhor - length(orig_20fullsample)))
)

write.csv(
  results_df,
  file.path(OUTPUT_DIR, "BDM_FAIR_robustness_probabilities.csv"),
  row.names = FALSE
)

message("Done. Outputs saved in: ", normalizePath(OUTPUT_DIR))
