polorth <- function(RR_SHOCK1, j) {
  r <- as.numeric(RR_SHOCK1)
  powers <- 0:j

  x <- outer(r, powers, `^`)
  x <- as.matrix(x)
  colnames(x) <- paste0("p", powers)

  list(x = x, D = 1)
}