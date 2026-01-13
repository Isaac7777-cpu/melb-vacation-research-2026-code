library(Matrix)
source(file = "lib/BFuns.R", verbose = TRUE)

# Negative log likelihood functions (up to multiplication 0.5 constant)
profile.nll <- function(
  eta,
  y,
  X = NULL,
  D,
  cov_model,
  resolution_matrix = NULL,
  resolution = NULL,
  nug = FALSE
) {
  if (is.null(resolution_matrix) && is.null(resolution)) {
    stop("You must supply iether `resolution_matrix` or at least `resolution`")
  }

  n <- length(y)

  if (is.null(resolution_matrix)) {
    resolution_matrix <- kronecker(
      Diagonal(n),
      Matrix(rep(1 / 9, 9), nrow = 1, sparse = TRUE)
    )
  }

  correlation_matrix <- cor.mat(D, eta, cov_model, nug)

  Q <- resolution_matrix %*% correlation_matrix %*% t(resolution_matrix)

  L <- t(chol(Q))
  log_det_cormat <- 2 * sum(log(diag(L)))

  white_1 <- solve(L, rep(1, n))
  white_y <- solve(L, y)

  # alpha_0 <- solve(t(white_1) %*% white_1, t(white_1) %*% white_y)
  # white_resids <- white_y - white_1 %*% alpha_0
  # No need to invert matrix to get the results
  alpha_0 <- (t(white_1) %*% white_y) / (t(white_1) %*% white_1)
  white_resids <- white_y - white_1 * as.numeric(alpha_0)

  mll2 <- log_det_cormat + n * log(mean(white_resids^2))
  return(mll2 / n)
}

MLE.fit <- function(
  y,
  D,
  cov_model,
  eta.ini,
  nug,
  opt,
  lo_bound,
  up_bound,
  negative_log_likelihood = profile.nll
) {
  m <- length(y)
  q <- length(eta.ini) + 1
  resolution_matrix <- kronecker(
    Diagonal(m),
    Matrix(rep(1 / 9, 9), nrow = 1, sparse = TRUE)
  )

  # Optimise using R base optimisation procedure--------------------------------
  if (opt == "LB") {
    prof.min <- optim(
      eta.ini,
      negative_log_likelihood,
      method = "L-BFGS-B",
      lower = lo_bound,
      upper = up_bound,
      y = y,
      cov_model = cov_model,
      D = D,
      nug = nug,
      resolution_matrix = resolution_matrix
    )
  } else if (opt == "NM") {
    prof.min <- optim(
      eta.ini,
      negative_log_likelihood,
      method = "Nelder-Mead",
      y = y,
      cov_model = cov_model,
      D = D,
      nug = nug,
      resolution_matrix = resolution_matrix
    )
  } else if (opt == "CG") {
    prof.min <- optim(
      eta.ini,
      negative_log_likelihood,
      method = "CG",
      y = y,
      cov_model = cov_model,
      D = D,
      nug = nug,
      resolution_matrix = resolution_matrix
    )
  }

  print("Finished optimising...")

  # Extract the results---------------------------------------------------------
  eta <- prof.min$par
  suc <- prof.min$convergence + 1
  nll <- prof.min$value

  # Calculate the estimated parameters------------------------------------------
  r_expanded <- cor.mat(D, eta, cov_model, nug)
  cormat <- resolution_matrix %*% r_expanded %*% t(resolution_matrix)
  L <- t(chol(cormat))

  white_1 <- solve(L, rep(1, m))
  white_y <- solve(L, y)
  alpha_0 <- (t(white_1) %*% white_y) / (t(white_1) %*% white_1)
  white_resids <- white_y - white_1 * as.numeric(alpha_0)

  # Calculate the variance of the estimated parameters--------------------------
  sill <- mean(white_resids^2)
  if (nug == TRUE) {
    nugget_effect <- sill * eta[q - 1]
    psill <- sill - nugget_effect
    theta_est <- c(eta[1:(q - 2)], psill, nugget_effect)
  } else if (nug == FALSE) {
    theta_est <- c(eta, sill)
  } else {
    stop("nug must be a boolean")
  }

  alpha_0_var <- m * sill

  # Save and return the results-------------------------------------------------
  # suc = 1 means convergence of optimization
  rlist <- list(
    theta = theta_est,
    suc,
    nll = nll,
    alpha0 = alpha_0,
    alpha_0_var = alpha_0_var,
    eta = eta
  )

  invisible(rlist)
}
