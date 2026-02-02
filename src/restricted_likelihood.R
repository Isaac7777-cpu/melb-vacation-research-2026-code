library(Matrix)
source(file = "lib/BFuns.R", verbose = TRUE)

profile.nll.reml.concentrated <- function(
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
    stop("You must supply either `resolution_matrix` or at least `resolution`")
  }

  n <- length(y)

  if (is.null(resolution_matrix)) {
    resolution_matrix <- kronecker(
      Diagonal(n),
      Matrix(rep(1 / resolution**2, resolution**2), nrow = 1, sparse = TRUE)
    )
  }

  correlation_matrix <- cor.mat(D, eta, cov_model, nug)
  Q <- resolution_matrix %*% correlation_matrix %*% t(resolution_matrix)
  L <- t(chol(Q))
  one <- rep(1, n)
  Linv_one <- forwardsolve(L, one)
  Linv_y <- forwardsolve(L, y)
  # alpha_0 <- (t(Linv_one) %*% Linv_y) / (t(Linv_one) %*% Linv_one)
  alpha_0 <- crossprod(Linv_one, Linv_y) / crossprod(Linv_one, Linv_one)
  Linv_resid <- Linv_y - as.numeric(alpha_0) * Linv_one
  # rss <- as.numeric(t(Linv_resid) %*% Linv_resid)
  rss <- as.numeric(crossprod(Linv_resid, Linv_resid))
  sigmasq <- rss / (n - 1)

  log_det_Q <- 2 * sum(log(diag(L)))

  nll2.reml <- (n - 1) +
    (n - 1) *
      log(2 * pi) +
    (n - 1) * log(sigmasq) +
    log_det_Q -
    log(n) +
    # as.numeric(log(t(Linv_one) %*% Linv_one))
    log(as.numeric(crossprod(Linv_one, Linv_one)))
  return(0.5 * nll2.reml)
}

REML.fit <- function(
  y,
  D,
  cov_model,
  eta.ini,
  nug,
  opt,
  lo_bound,
  up_bound,
  negative_log_likelihood = profile.nll.reml.concentrated,
  resolution,
  verbose = FALSE
) {
  m <- length(y)
  q <- length(eta.ini) + 1
  resolution_matrix <- kronecker(
    Diagonal(m),
    Matrix(rep(1 / resolution^2, resolution^2), nrow = 1, sparse = TRUE)
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

  if (verbose) {
    print("Finished optimising...")
  }

  # Extract the results---------------------------------------------------------
  eta <- prof.min$par
  suc <- prof.min$convergence + 1
  nll <- prof.min$value

  # Calculate the estimated parameters------------------------------------------
  correlation_matrix <- cor.mat(D, eta, cov_model, nug)
  Q <- resolution_matrix %*% correlation_matrix %*% t(resolution_matrix)
  L <- t(chol(Q))
  one <- rep(1, m)
  Linv_one <- forwardsolve(L, one)
  Linv_y <- forwardsolve(L, y)
  alpha_0 <- (t(Linv_one) %*% Linv_y) / (t(Linv_one) %*% Linv_one)
  Linv_resid <- Linv_y - as.numeric(alpha_0) * Linv_one
  rss <- as.numeric(t(Linv_resid) %*% Linv_resid)
  sill <- rss / (m - 1)

  # Calculate the variance of the estimated parameters--------------------------
  if (nug == TRUE) {
    nugget_effect <- sill * eta[q - 1]
    psill <- sill - nugget_effect
    theta_est <- c(eta[1:(q - 2)], psill, nugget_effect)
  } else if (nug == FALSE) {
    theta_est <- c(eta, sill)
  } else {
    stop("nug must be a boolean")
  }

  # alpha_0_var <- m * sill
  alpha_0_var <- sill / as.numeric(crossprod(Linv_one, Linv_one))

  # Save and return the results-------------------------------------------------
  # suc = 1 means convergence of optimization
  rlist <- list(
    theta = theta_est,
    "converge (1 = yes)" = suc,
    nll = nll,
    alpha0 = alpha_0,
    alpha_0_var = alpha_0_var,
    eta = eta
  )

  invisible(rlist)
}

profile.nll.reml.concentrated.memeff <- function(
  eta,
  y,
  X = NULL,
  coord,
  cov_model,
  resolution,
  nug = FALSE
) {
  n <- length(y)

  Q <- build_Q_blockavg(
    coord = coord,
    resolution = resolution,
    eta = eta,
    cov_model = cov_model,
    nug = nug
  )
  L <- t(chol(Q))
  one <- rep(1, n)
  Linv_one <- forwardsolve(L, one)
  Linv_y <- forwardsolve(L, y)
  # alpha_0 <- (t(Linv_one) %*% Linv_y) / (t(Linv_one) %*% Linv_one)
  alpha_0 <- crossprod(Linv_one, Linv_y) / crossprod(Linv_one, Linv_one)
  Linv_resid <- Linv_y - as.numeric(alpha_0) * Linv_one
  # rss <- as.numeric(t(Linv_resid) %*% Linv_resid)
  rss <- as.numeric(crossprod(Linv_resid, Linv_resid))
  sigmasq <- rss / (n - 1)

  log_det_Q <- 2 * sum(log(diag(L)))

  nll2.reml <- (n - 1) +
    (n - 1) *
      log(2 * pi) +
    (n - 1) * log(sigmasq) +
    log_det_Q -
    log(n) +
    # as.numeric(log(t(Linv_one) %*% Linv_one))
    log(as.numeric(crossprod(Linv_one, Linv_one)))
  return(0.5 * nll2.reml)
}

REML.fit.memeff <- function(
  y,
  coord,
  cov_model,
  eta.ini,
  nug,
  opt,
  lo_bound,
  up_bound,
  negative_log_likelihood = profile.nll.reml.concentrated.memeff,
  resolution,
  verbose = FALSE
) {
  m <- length(y)
  q <- length(eta.ini) + 1

  # Optimise using R base optimisation procedure--------------------------------
  if (opt == "LB") {
    prof.min <- optim(
      eta.ini,
      negative_log_likelihood,
      method = "L-BFGS-B",
      lower = lo_bound,
      upper = up_bound,
      y = y,
      coord = coord,
      cov_model = cov_model,
      nug = nug,
      resolution = resolution,
    )
  } else if (opt == "NM") {
    prof.min <- optim(
      eta.ini,
      negative_log_likelihood,
      method = "Nelder-Mead",
      y = y,
      coord = coord,
      cov_model = cov_model,
      nug = nug,
      resolution = resolution
    )
  } else if (opt == "CG") {
    prof.min <- optim(
      eta.ini,
      negative_log_likelihood,
      method = "CG",
      y = y,
      coord = coord,
      cov_model = cov_model,
      nug = nug,
      resolution = resolution
    )
  }

  if (verbose) {
    print("Finished optimising...")
  }

  # Extract the results---------------------------------------------------------
  eta <- prof.min$par
  suc <- prof.min$convergence + 1
  nll <- prof.min$value

  # Calculate the estimated parameters------------------------------------------
  Q <- build_Q_blockavg(
    coord = coord,
    resolution = resolution,
    eta = eta,
    cov_model = cov_model,
    nug = nug
  )
  L <- t(chol(Q))
  one <- rep(1, m)
  Linv_one <- forwardsolve(L, one)
  Linv_y <- forwardsolve(L, y)
  alpha_0 <- (t(Linv_one) %*% Linv_y) / (t(Linv_one) %*% Linv_one)
  Linv_resid <- Linv_y - as.numeric(alpha_0) * Linv_one
  rss <- as.numeric(t(Linv_resid) %*% Linv_resid)
  sill <- rss / (m - 1)

  # Calculate the variance of the estimated parameters--------------------------
  if (nug == TRUE) {
    nugget_effect <- sill * eta[q - 1]
    psill <- sill - nugget_effect
    theta_est <- c(eta[1:(q - 2)], psill, nugget_effect)
  } else if (nug == FALSE) {
    theta_est <- c(eta, sill)
  } else {
    stop("nug must be a boolean")
  }

  # alpha_0_var <- m * sill
  alpha_0_var <- sill / as.numeric(crossprod(Linv_one, Linv_one))

  # Save and return the results-------------------------------------------------
  # suc = 1 means convergence of optimization
  rlist <- list(
    theta = theta_est,
    "converge (1 = yes)" = suc,
    nll = nll,
    alpha0 = alpha_0,
    alpha_0_var = alpha_0_var,
    eta = eta
  )

  invisible(rlist)
}
