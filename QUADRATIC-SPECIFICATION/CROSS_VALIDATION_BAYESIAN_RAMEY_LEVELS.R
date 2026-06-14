source("polorth.R")
source("nwest_balanced_panel_bayesian.R")
source("generate_unobserved_shock_ramey_levels.R")
source("bayes_utils.R")

get_mat_var <- function(mat, name) {
  if (!is.null(mat[[name]])) return(mat[[name]])
  dotted <- gsub("_", ".", name)
  if (!is.null(mat[[dotted]])) return(mat[[dotted]])
  lower_names <- tolower(gsub("[._]", "", names(mat)))
  target <- tolower(gsub("[._]", "", name))
  idx <- match(target, lower_names)
  if (!is.na(idx)) return(mat[[idx]])
  NULL
}

CROSS_VALIDATION_BAYESIAN_RAMEY_LEVELS <- function(u, mon, oil, fis, pol,
                                                    cv_ts, cv_ts_block,
                                                    aic_shock, aic, mc,
                                                    IP, RR_SHOCK, unc,
                                                    lag_ramey, nlag,
                                                    data_file = "FISCAL_DATA_SHOCK.mat",
                                                    rep = 2000L,
                                                    seed = NULL) {
  # R translation focused on the levels/fiscal branch used by FINAL_OUTPUT.m.
  if (!is.null(seed)) set.seed(seed)

  if (aic != 1) {
    stop("This translated function currently implements the aic == 1 branch used by FINAL_OUTPUT.m.")
  }
  if (!(fis == 1 && mon == 0 && oil == 0 && unc == 0)) {
    stop("This translated function implements the fiscal branch used by FINAL_OUTPUT.m.")
  }
  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    stop("Package 'R.matlab' is required to read .mat files. Install with install.packages('R.matlab').")
  }

  ee <- 0L
  h <- 1:20

  if (mc == 0) {
    mat <- R.matlab::readMat(data_file)
    IP <- get_mat_var(mat, "IP")
    RR_SHOCK <- get_mat_var(mat, "RR_SHOCK")
    if (is.null(IP) || is.null(RR_SHOCK)) stop("Could not read IP and RR_SHOCK from data_file.")
    if (lag_ramey == 0) {
      IP1 <- IP[(nlag + 1):nrow(IP), , drop = FALSE]
    } else {
      IP1 <- IP[(2 * nlag + 1):nrow(IP), , drop = FALSE]
    }
  } else {
    IP1 <- IP
  }

  H <- max(h)
  lead_coeff_y <- array(0, dim = c(H, pol, rep))
  lead_coeff_g <- array(0, dim = c(H, pol, rep))
  lead_coeff_y1 <- array(0, dim = c(H, pol, rep))
  lead_coeff_g1 <- array(0, dim = c(H, pol, rep))

  e_y <- matrix(vector("list", H * rep), nrow = H, ncol = rep)
  e_g <- matrix(vector("list", H * rep), nrow = H, ncol = rep)
  reg_y <- array(vector("list", rep * H * pol), dim = c(rep, H, pol))
  reg_g <- array(vector("list", rep * H * pol), dim = c(rep, H, pol))

  xx <- polorth(1, pol)$x
  xx1 <- polorth(-1, pol)$x
  xx2 <- matrix(c(1, 0, 0), nrow = 1)
  if (ncol(xx2) < pol + 1) xx2 <- cbind(xx2, matrix(0, 1, pol + 1 - ncol(xx2)))

  for (t in seq_len(rep)) {
    if (lag_ramey == 0) {
      gen <- GENERATE_UNOBSERVED_SHOCK_RAMEY_LEVELS(1, 0, 0, 0, nlag, data_file = data_file)
      RR_SHOCK1 <- gen$resid_final
      x <- polorth(RR_SHOCK1, pol)$x
      s <- 0L
    } else {
      gen <- GENERATE_UNOBSERVED_SHOCK_RAMEY_LEVELS(1, 0, 0, 0, nlag, data_file = data_file)
      RR_SHOCK1 <- as.numeric(gen$resid_final)
      ip2 <- gen$ip2
      xpowers <- sapply(1:pol, function(p) RR_SHOCK1^p)
      xpowers <- as.matrix(xpowers)

      # MATLAB: x1 = lagmatrix(x, 0:nlag-1); x1 = x1(nlag+1:end,:)
xlag_list <- vector("list", nlag)

for (lag in 0:(nlag - 1)) {
  if (lag == 0) {
    xlag_list[[lag + 1]] <- xpowers
  } else {
    xlag_list[[lag + 1]] <- rbind(
      matrix(NA_real_, nrow = lag, ncol = ncol(xpowers)),
      xpowers[1:(nrow(xpowers) - lag), , drop = FALSE]
    )
  }
}

    x1 <- do.call(cbind, xlag_list)
    x1 <- x1[(nlag + 1):nrow(x1), , drop = FALSE]

    ip2_trim <- ip2[(nlag + 1):nrow(ip2), , drop = FALSE]

    common_n <- min(nrow(x1), nrow(ip2_trim))
    x1 <- x1[seq_len(common_n), , drop = FALSE]
    ip2_trim <- ip2_trim[seq_len(common_n), , drop = FALSE]

    x <- cbind(1, x1, ip2_trim)
    s <- ncol(cbind(x1[, (pol + 1):ncol(x1), drop = FALSE], ip2_trim))
    }

    for (i in h) {
      Tn <- nrow(IP1[(i):nrow(IP1), 1, drop = FALSE])
      y <- rbind(matrix(IP1[i:nrow(IP1), 1], ncol = 1),
                 matrix(IP1[i:nrow(IP1), 2], ncol = 1))

      for (j in seq_len(pol)) {
        if (is.null(e_y[[i, t]])) e_y[[i, t]] <- matrix(NA_real_, Tn, pol)
        if (is.null(e_g[[i, t]])) e_g[[i, t]] <- matrix(NA_real_, Tn, pol)
        rows_x <- 1:(nrow(x) - i + 1 - ee)
        if (lag_ramey == 0) {
          cols <- 1:(j + 1)
          xb <- x[rows_x, cols, drop = FALSE]
          Xbig <- rbind(cbind(xb, matrix(0, nrow(xb), ncol(xb))),
                        cbind(matrix(0, nrow(xb), ncol(xb)), xb))
        } else {
          cols <- c(1:(j + 1), (ncol(x) - s + 1):ncol(x))
          xb <- x[rows_x, cols, drop = FALSE]
          Xbig <- rbind(cbind(xb, matrix(0, nrow(xb), ncol(xb))),
                        cbind(matrix(0, nrow(xb), ncol(xb)), xb))
        }

        results <- nwest_balanced_panel_bayesian(y, Xbig, i, Tn, 2)
        b <- results$beta

        vneg <- matrix(xx1[1, 1:(j + 1)], ncol = 1)
        vzero <- matrix(xx2[1, 1:(j + 1)], ncol = 1)
        vpos <- matrix(xx[1, 1:(j + 1)], ncol = 1)

        by <- b[1:(j + 1), , drop = FALSE]
        lead_coeff_y[i, j, t] <- -as.numeric(t(by) %*% vneg - t(by) %*% vzero)
        lead_coeff_y1[i, j, t] <- as.numeric(t(by) %*% vpos - t(by) %*% vzero)
        e_y[[i, t]][, j] <- (IP1[i:nrow(IP1), 1] - x[rows_x, 1:(j + 1), drop = FALSE] %*% by) / i
        reg_y[[t, i, j]] <- by

        if (lag_ramey == 0) {
          bg <- b[(j + 2):nrow(b), , drop = FALSE]
        } else {
          bg <- b[(j + 2 + s):(nrow(b) - s), , drop = FALSE]
        }
        lead_coeff_g[i, j, t] <- -as.numeric(t(bg) %*% vneg - t(bg) %*% vzero)
        lead_coeff_g1[i, j, t] <- as.numeric(t(bg) %*% vpos - t(bg) %*% vzero)
        e_g[[i, t]][, j] <- (IP1[i:nrow(IP1), 2] - x[rows_x, 1:(j + 1), drop = FALSE] %*% bg) / i
        reg_g[[t, i, j]] <- bg
      }
    }
  }

  list(
    LEAD_COEFF_Y = lead_coeff_y,
    LEAD_COEFF_G = lead_coeff_g,
    residual_y = e_y[1, , drop = FALSE],
    residual_g = e_g[1, , drop = FALSE],
    reg_y = reg_y,
    reg_g = reg_g,
    LEAD_COEFF_Y1 = lead_coeff_y1,
    LEAD_COEFF_G1 = lead_coeff_g1
  )
}
