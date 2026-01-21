#' Calculate the variance
cov.exponential <- function(dist, sigmasq, phi, tausq = NA) {
  cov.mat <- sigmasq * exp(-dist / phi)
  if (!is.na(tausq)) {
    diag(cov.mat) <- diag(cov.mat) + tausq
  }
  cov.mat
}
