# Translation of FINAL_OUTPUT.m for the baseline levels fiscal specification.
# Required packages: install.packages(c("R.matlab", "MASS"))

source("bayes_utils.R")
source("CROSS_VALIDATION_BAYESIAN_RAMEY_LEVELS.R")

run_final_output <- function(data_file = "FISCAL_DATA_SHOCK.mat",
                             out_dir = ".",
                             rep = 2000L,
                             seed = NULL) {
  if (!dir.exists("QUAD_R_OUTPUT")) {
  dir.create("QUAD_R_OUTPUT", recursive = TRUE)
}
  cv_ts <- 0; cv_ts_block <- 0; aic_shock <- 0; aic <- 1
  recursive <- 0
  levels <- 1
  state <- 0
  lag_ramey <- 1
  nlag <- 5
  mc <- 0
  a <- 0.025
  pol <- if (aic_shock == 1) 1 else 2
  RR_SHOCK <- 0
  IP <- 0

  if (state != 0 || levels != 1 || recursive != 0) {
    stop("This translated script implements the state=0, levels=1, recursive=0 branch in FINAL_OUTPUT.m.")
  }

  CV_FIS <- CROSS_VALIDATION_BAYESIAN_RAMEY_LEVELS(
    u = 1, mon = 0, oil = 0, fis = 1, pol = pol,
    cv_ts = cv_ts, cv_ts_block = cv_ts_block,
    aic_shock = aic_shock, aic = aic, mc = mc,
    IP = IP, RR_SHOCK = RR_SHOCK, unc = 0,
    lag_ramey = lag_ramey, nlag = nlag,
    data_file = data_file, rep = rep, seed = seed
  )

  L <- dim(CV_FIS$LEAD_COEFF_G)[2]
  L1 <- dim(CV_FIS$LEAD_COEFF_G)[3]
  L4 <- dim(CV_FIS$reg_g)[2]

  reg_y_list <- vector("list", L)
  reg_g_list <- vector("list", L)
  for (i in seq_len(L)) {
    # Each entry is coefficient x horizon x replication.
    first_y <- CV_FIS$reg_y[[1, 1, i]]
    first_g <- CV_FIS$reg_g[[1, 1, i]]
    reg_y_arr <- array(NA_real_, dim = c(nrow(first_y), L4, L1))
    reg_g_arr <- array(NA_real_, dim = c(nrow(first_g), L4, L1))
    for (t in seq_len(L1)) {
      for (h in seq_len(L4)) {
        reg_y_arr[, h, t] <- as.numeric(CV_FIS$reg_y[[t, h, i]])
        reg_g_arr[, h, t] <- as.numeric(CV_FIS$reg_g[[t, h, i]])
      }
    }
    reg_y_list[[i]] <- reg_y_arr
    reg_g_list[[i]] <- reg_g_arr
  }

  if (L1 > 999) {
    probs <- c(a, 0.5, 1 - a)
    irf_y_pos <- col_quantile(CV_FIS$LEAD_COEFF_Y[, 2, , drop = FALSE], probs)
    irf_y_neg <- col_quantile(CV_FIS$LEAD_COEFF_Y1[, 2, , drop = FALSE], probs)
    irf_y_asym <- col_quantile(CV_FIS$LEAD_COEFF_Y1[, 2, , drop = FALSE] - CV_FIS$LEAD_COEFF_Y[, 2, , drop = FALSE], probs)
    irf_y_lin <- col_quantile(CV_FIS$LEAD_COEFF_Y1[, 1, , drop = FALSE], probs)

    irf_g_pos <- col_quantile(CV_FIS$LEAD_COEFF_G[, 2, , drop = FALSE], probs)
    irf_g_neg <- col_quantile(CV_FIS$LEAD_COEFF_G1[, 2, , drop = FALSE], probs)
    irf_g_asym <- col_quantile(CV_FIS$LEAD_COEFF_G1[, 2, , drop = FALSE] - CV_FIS$LEAD_COEFF_G[, 2, , drop = FALSE], probs)
    irf_g_lin <- col_quantile(CV_FIS$LEAD_COEFF_G1[, 1, , drop = FALSE], probs)

    cy_pos <- cumsum_dim1(CV_FIS$LEAD_COEFF_Y[, 2, , drop = FALSE])
    cg_pos <- cumsum_dim1(CV_FIS$LEAD_COEFF_G[, 2, , drop = FALSE])
    cy_neg <- cumsum_dim1(CV_FIS$LEAD_COEFF_Y1[, 2, , drop = FALSE])
    cg_neg <- cumsum_dim1(CV_FIS$LEAD_COEFF_G1[, 2, , drop = FALSE])
    cy_lin <- cumsum_dim1(CV_FIS$LEAD_COEFF_Y[, 1, , drop = FALSE])
    cg_lin <- cumsum_dim1(CV_FIS$LEAD_COEFF_G[, 1, , drop = FALSE])

    irf_m_pos <- col_quantile(cy_pos / cg_pos, probs)
    irf_m_neg <- col_quantile(cy_neg / cg_neg, probs)
    irf_m_asym <- col_quantile((cy_neg / cg_neg) - (cy_pos / cg_pos), probs)
    irf_m_lin <- col_quantile(cy_lin / cg_lin, probs)

    reg_y <- apply(reg_y_list[[2]], c(1, 2), stats::median, na.rm = TRUE)
    reg_g <- apply(reg_g_list[[2]], c(1, 2), stats::median, na.rm = TRUE)
    residual_y <- CV_FIS$residual_y[, 2, drop = FALSE]
    residual_g <- CV_FIS$residual_g[, 2, drop = FALSE]
    num_g <- sum(irf_g_asym[, , 1] > 0 | irf_g_asym[, , 3] < 0)
    num_y <- sum(irf_y_asym[, , 1] > 0 | irf_y_asym[, , 3] < 0)

    save_file <- file.path(out_dir, sprintf("FIS_MC_INPUT_LEVELS_%d.RData", lag_ramey))
    save(reg_y, reg_g, residual_g, residual_y,
         irf_y_neg, irf_y_pos, irf_y_asym,
         irf_g_neg, irf_g_pos, irf_g_asym,
         irf_m_neg, irf_m_pos, irf_m_asym, irf_m_lin,
         num_g, num_y, file = save_file)
  } else {
    stop("The original MATLAB plotting/saving branch expects at least 1000 posterior replications.")
  }

  plot_irfs(out_dir, lag_ramey,
            irf_y_neg, irf_y_pos, irf_y_asym,
            irf_g_neg, irf_g_pos, irf_g_asym,
            irf_m_neg, irf_m_pos, irf_m_asym,
            irf_y_lin, irf_g_lin, irf_m_lin)

  invisible(list(CV_FIS = CV_FIS,
                 irfs = list(irf_y_neg = irf_y_neg, irf_y_pos = irf_y_pos,
                             irf_y_asym = irf_y_asym, irf_g_neg = irf_g_neg,
                             irf_g_pos = irf_g_pos, irf_g_asym = irf_g_asym,
                             irf_m_neg = irf_m_neg, irf_m_pos = irf_m_pos,
                             irf_m_asym = irf_m_asym, irf_y_lin = irf_y_lin,
                             irf_g_lin = irf_g_lin, irf_m_lin = irf_m_lin)))
}

plot_band <- function(x, arr, scale = 1, main = "", ylab = "") {
  plot(x, scale * arr[, 1, 2], type = "l", lwd = 3, lty = 1,
       xlab = "Horizon", ylab = ylab, main = main)
  lines(x, scale * arr[, 1, 1], lwd = 3, lty = 2)
  lines(x, scale * arr[, 1, 3], lwd = 3, lty = 2)
}

plot_irfs <- function(out_dir, lag_ramey,
                      irf_y_neg, irf_y_pos, irf_y_asym,
                      irf_g_neg, irf_g_pos, irf_g_asym,
                      irf_m_neg, irf_m_pos, irf_m_asym,
                      irf_y_lin, irf_g_lin, irf_m_lin) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  x <- 1:20
  w <- 100

  pdf(file.path(out_dir, sprintf("FISCAL_NEWS_FIGURE_NONLINEAR_LEVELS1%d.pdf", lag_ramey)), width = 12, height = 9)
  par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))
  plot_band(x, irf_y_neg, w, "GDP: Positive Shock", "Percentage")
  plot_band(x, irf_y_pos, w, "GDP: Negative Shock", "Percentage")
  plot_band(x, irf_y_asym, w, "GDP: Response Asymmetry", "Percentage Points")
  plot_band(x, irf_g_neg, w, "Government Spending: Positive Shock", "Percentage")
  plot_band(x, irf_g_pos, w, "Government Spending: Negative Shock", "Percentage")
  plot_band(x, irf_g_asym, w, "Government Spending: Response Asymmetry", "Percentage Points")
  plot_band(x, irf_m_neg, 1, "Fiscal Multiplier: Positive Shock", "Fiscal Multiplier Units")
  plot_band(x, irf_m_pos, 1, "Fiscal Multiplier: Negative Shock", "Fiscal Multiplier Units")
  plot_band(x, irf_m_asym, 1, "Fiscal Multiplier: Asymmetry", "Fiscal Multiplier Units")
  dev.off()

  pdf(file.path(out_dir, sprintf("FISCAL_NEWS_FIGURE_LINEAR_LEVELS1%d.pdf", lag_ramey)), width = 12, height = 8)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  plot_band(x, irf_y_lin, w, "GDP: Linear Specification", "Percentage")
  plot_band(x, irf_g_lin, w, "Government Spending: Linear Specification", "Percentage")
  plot_band(x, irf_m_lin, 1, "Fiscal Multiplier: Linear Specification", "Fiscal Multiplier Units")
  plot.new()
  dev.off()
}

run_final_output(data_file = "FISCAL_DATA_SHOCK.mat", out_dir = "QUAD_R_OUTPUT", rep = 2000, seed = 123)
# Example:
# setwd("/path/to/translation_folder")
# run_final_output(data_file = "../FISCAL_DATA_SHOCK.mat", out_dir = "R_OUTPUT", rep = 2000, seed = 123)
