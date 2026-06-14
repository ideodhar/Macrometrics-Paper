source("bayes_utils.R")

nwest_balanced_panel_bayesian <- function(yyy, x, nlag, T, N) {
  # Translation of nwest_balanced_panel_bayesian.m
  yyy <- matrix(as.numeric(yyy), ncol = 1)
  x <- as.matrix(x)
  nobs_total <- nrow(x)
  nvar <- ncol(x)
  nobs <- nobs_total / N

  xpxi <- safe_solve(crossprod(x))
  beta <- xpxi %*% crossprod(x, yyy)
  yhat <- x %*% beta
  resid <- yyy - yhat
  sigu <- as.numeric(crossprod(resid))

  # Newey-West correction, following the MATLAB code literally.
  emat <- matrix(rep(as.numeric(resid), each = nvar), nrow = nvar, byrow = FALSE)
  hhat <- emat * t(x)
  G <- matrix(0, nvar, nvar)
  uu <- matrix(NA_real_, nvar, nobs)

  for (i in seq_len(nvar)) {
    u <- array(hhat[i, ], dim = c(T, N))
    u1 <- t(u)
    uu[i, ] <- colSums(u1)
  }

  for (a in 0:nlag) {
    weight <- (nlag + 1 - a) / (nlag + 1)
    za <- uu[, (a + 1):nobs, drop = FALSE] %*% t(uu[, 1:(nobs - a), drop = FALSE])
    ga <- if (a == 0) za else za + t(za)
    G <- G + weight * ga
  }

  V <- as.numeric(sum(diag(xpxi %*% G)))
  SIG <- iwishrnd_scalar(sigu + V, sum(resid != 0))
  beta_draw <- mvnrnd_vec(as.numeric(beta), SIG * xpxi)

  list(
    meth = "nwest",
    yyy = yyy,
    nobs = nobs,
    nvar = nvar,
    beta = matrix(beta_draw, ncol = 1),
    yhat = yhat,
    resid = resid
  )
}
