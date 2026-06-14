iwishrnd_scalar <- function(scale, df) {
  # Scalar inverse-Wishart draw. MATLAB iwishrnd(S, df) when S is scalar.
  as.numeric(scale) / stats::rchisq(1, df = df)
}

mvnrnd_vec <- function(mu, Sigma) {
  mu <- as.numeric(mu)
  Sigma <- as.matrix(Sigma)
  if (length(mu) == 1L) {
    return(stats::rnorm(1, mean = mu, sd = sqrt(as.numeric(Sigma))))
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for multivariate normal draws.")
  }
  as.numeric(MASS::mvrnorm(n = 1, mu = mu, Sigma = Sigma))
}

safe_solve <- function(A, b = NULL) {
  A <- as.matrix(A)
  if (is.null(b)) {
    return(solve(A))
  }
  solve(A, b)
}

col_quantile <- function(arr, probs) {
  # MATLAB: quantile(array, probs, 3) for array with dimensions H x 1 x Rep.
  d <- dim(arr)
  if (length(d) != 3L) stop("Expected a 3D array.")
  out <- array(NA_real_, dim = c(d[1], d[2], length(probs)))
  for (i in seq_len(d[1])) {
    for (j in seq_len(d[2])) {
      out[i, j, ] <- stats::quantile(arr[i, j, ], probs = probs, na.rm = TRUE, names = FALSE, type = 7)
    }
  }
  out
}

cumsum_dim1 <- function(arr) {
  d <- dim(arr)
  out <- arr
  for (j in seq_len(d[2])) {
    for (k in seq_len(d[3])) {
      out[, j, k] <- cumsum(arr[, j, k])
    }
  }
  out
}
