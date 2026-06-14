# Figure6_translated.R
# Translation of Figure6.m from MATLAB to R
# Requires BRZDAT.xlsx in the working directory, with the US data on sheet 2.

rm(list = ls())
graphics.off()

# Packages
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required. Install it with: install.packages('readxl')")
}

# -----------------------------------------------------------------------------
# DATA
# -----------------------------------------------------------------------------
vdata <- as.matrix(readxl::read_excel("BRZDAT.xlsx", sheet = 2, col_names = FALSE))
storage.mode(vdata) <- "numeric"

time <- vdata[, 1]

# -----------------------------------------------------------------------------
# IMPORTANT INPUT: sample window
# -----------------------------------------------------------------------------
t1 <- which(time == 1890)
t2 <- which(time == 2015.75)

if (length(t1) == 0 || length(t2) == 0) {
  stop("Could not find sample endpoints 1890 and/or 2015.75 in column 1 of BRZDAT.xlsx.")
}

vdata <- vdata[t1:t2, , drop = FALSE]

q         <- vdata[, 1]
ngov      <- vdata[, 2]
ngdp      <- vdata[, 3]
pgdp      <- vdata[, 4]
totpop    <- vdata[, 5]
recession <- vdata[, 6]
unemp     <- vdata[, 7]
pdvmil    <- vdata[, 8]
realgdp   <- vdata[, 9]
ntax      <- vdata[, 10]
potgdp    <- vdata[, 11]

# -----------------------------------------------------------------------------
# DATA TRANSFORMATION
# -----------------------------------------------------------------------------
# MATLAB behavior: pdvmily(1) is implicitly zero because the loop starts at t = 2.
pdvmily <- numeric(length(pdvmil))
for (t in 2:length(pdvmil)) {
  pdvmily[t] <- pdvmil[t] / (pgdp[t - 1] * potgdp[t - 1])
}

rgdp     <- realgdp / totpop
rgov     <- ngov / totpop / pgdp
lrgdp    <- log(realgdp / totpop)
lrgov    <- log(ngov / totpop / pgdp)
lrtax    <- log(ntax / totpop / pgdp)
rgdp_pot <- realgdp / potgdp
rgov_pot <- ngov / potgdp / pgdp

data <- cbind(pdvmily * 100, rgdp_pot, rgov_pot, unemp)

varnmz <- c(
  "Military news (% of GDP)",
  "Normalized GDP",
  "Normalized govt. spending",
  "Unemployment rate"
)

# -----------------------------------------------------------------------------
# FIGURE 6
# -----------------------------------------------------------------------------
# Equivalent of MATLAB vline() calls using abline(v = ...).
dir.create("R_OUTPUT", showWarnings = FALSE)

pdf("R_OUTPUT/figure_6.pdf", width = 10, height = 7)
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par), add = TRUE)

par(mfrow = c(2, 2), lwd = 1.5)

for (i in 1:4) {
  y <- data[, i]

  plot(
    q, y,
    type = "l",
    col = "blue",
    lwd = 1,
    xlim = c(q[1], q[length(q)]),
    ylim = c(min(y, na.rm = TRUE) - 0.05, max(y, na.rm = TRUE) + 0.05),
    xlab = "",
    ylab = "",
    main = varnmz[i]
  )

  abline(v = 1901.25, col = "red", lty = 2)
  abline(v = 2015.75, col = "red", lty = 2)
  abline(v = 1890,    col = "black", lty = 2)
  abline(v = 2014.75, col = "black", lty = 2)
}

dev.off()
