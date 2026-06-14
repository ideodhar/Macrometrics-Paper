source("bayes_utils.R")

GENERATE_UNOBSERVED_SHOCK_RAMEY_LEVELS <- function(tfp, cs, unc, mon, nlag,
                                                    data_file = "FISCAL_DATA_SHOCK.mat") {
  # Translation of GENERATE_UNOBSERVED_SHOCK_RAMEY_LEVELS.m for the fiscal/Ramey case.
  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    stop("Package 'R.matlab' is required to read .mat files. Install with install.packages('R.matlab').")
  }

  if (tfp != 1 || cs != 0 || unc != 0 || mon != 0) {
    stop("This R translation implements the fiscal/Ramey levels branch used by FINAL_OUTPUT.m: tfp=1, cs=0, unc=0, mon=0.")
  }

  mat <- R.matlab::readMat(data_file)
  IP <- mat$IP
  RR_SHOCK <- mat$RR.SH0CK
  if (is.null(RR_SHOCK)) RR_SHOCK <- mat$RR_SHOCK
  if (is.null(RR_SHOCK)) RR_SHOCK <- mat$RR.SH0CK
  if (is.null(RR_SHOCK)) {
    # R.matlab converts underscores to dots by default.
    RR_SHOCK <- mat$RR.SHOCK
  }
  if (is.null(IP) || is.null(RR_SHOCK)) stop("Could not find IP and RR_SHOCK in FISCAL_DATA_SHOCK.mat.")

  IP <- IP[2:nrow(IP), , drop = FALSE]
  F <- as.matrix(RR_SHOCK)
  ip <- IP[, 1:2, drop = FALSE]

  # MATLAB loops j = 1:(nlag-1), building lagged shock and lagged IP block.
  n <- nrow(F) - nlag + 1
  F1 <- matrix(NA_real_, nrow = n, ncol = nlag - 1)
  ip_blocks <- vector("list", nlag - 1)
  for (j in 1:(nlag - 1)) {
    rows <- (nlag - j):(nrow(F) - j)
    F1[, j] <- F[rows, 1]
    ip_blocks[[j]] <- ip[rows, , drop = FALSE]
  }
  ip2 <- do.call(cbind, ip_blocks)

  x <- cbind(F1, ip2, 1)
  yyy <- F[nlag:nrow(F), 1]
  z <- cbind(x, yyy)
  z <- z[stats::complete.cases(z), , drop = FALSE]
  yyy <- matrix(z[, ncol(z)], ncol = 1)
  x1 <- as.matrix(z[, -ncol(z), drop = FALSE])

  xpxi <- safe_solve(crossprod(x1))
  beta <- xpxi %*% crossprod(x1, yyy)
  yhat <- x1 %*% beta
  resid <- yyy - yhat
  sigu <- as.numeric(crossprod(resid))

  SIG <- iwishrnd_scalar(sigu, sum(resid != 0))
  beta <- matrix(mvnrnd_vec(as.numeric(beta), SIG * xpxi), ncol = 1)
  yhat <- x1 %*% beta
  resid <- yyy - yhat
  resid_final <- as.numeric(resid) / sqrt(SIG)

  list(resid_final = matrix(resid_final, ncol = 1), ip2 = ip2)
}
