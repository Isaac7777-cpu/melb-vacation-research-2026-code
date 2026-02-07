# This function corresponds to MLE_fit and MLE_fitE in Python

MLE.fit <- function(y, X, D, cov_model, eta.ini, nug, opt, lo_bound, up_bound) {
  # Preliminaries:
  n <- length(y)
  q <- length(eta.ini) + 1
  # minimze negative concentrated log-likelihood function

  if (opt == "LB") {
    prof.min <- optim(
      eta.ini,
      profile.nll,
      method = "L-BFGS-B",
      lower = lo_bound,
      upper = up_bound,
      y = y,
      X = X,
      cov_model = cov_model,
      D = D,
      nug = nug
    )
  } else if (opt == "NM") {
    prof.min <- optim(
      eta.ini,
      profile.nll,
      method = "Nelder-Mead",
      y = y,
      X = X,
      cov_model = cov_model,
      D = D,
      nug = nug
    )
  } else if (opt == "CG") {
    prof.min <- optim(
      eta.ini,
      profile.nll,
      method = "CG",
      y = y,
      X = X,
      cov_model = cov_model,
      D = D,
      nug = nug
    )
  }

  # Return theta
  eta <- prof.min$par
  suc <- prof.min$convergence + 1
  nll <- prof.min$value

  cormat <- cor.mat(D, eta, cov_model, nug)
  L <- t(chol(cormat))

  if (is.null(X)) {
    white_resids <- solve(L, y)
  } else {
    white_X <- solve(L, X)
    white_y <- solve(L, y)

    beta_est <- solve(t(white_X) %*% white_X, t(white_X) %*% white_y)
    white_resids <- white_y - white_X %*% beta_est
  }

  sill = mean(white_resids^2)
  if (nug == TRUE) {
    nugget_effect <- sill * eta[q - 1]
    psill <- sill - nugget_effect
    theta_est <- c(eta[1:(q - 2)], psill, nugget_effect)
  } else if (nug == FALSE) {
    theta_est <- c(eta, sill)
  }

  if (!is.null(X)) {
    beta_var <- solve(t(white_X) %*% white_X) * sill
  }

  # suc = 1 means convergence of optimization
  if (is.null(X)) {
    rlist <- list(theta = theta_est, suc, nll = nll, eta = eta, sill = sill)
  } else {
    rlist <- list(
      theta = theta_est,
      suc,
      nll = nll,
      beta = beta_est,
      beta_var = beta_var,
      eta = eta
    )
  }

  invisible(rlist)
}


# This function correpsonds to profile_nll and profile_nllE in Python.
profile.nll <- function(eta, y, X = NULL, D, cov_model, nug) {
  n <- length(y)
  cormat <- cor.mat(D, eta, cov_model, nug)

  L <- t(chol(cormat))
  log_det_cormat = 2 * sum(log(diag(L)))

  if (is.null(X)) {
    white_resids = solve(L, y)
  } else {
    white_X <- solve(L, X)
    white_y <- solve(L, y)

    beta <- solve(t(white_X) %*% white_X, t(white_X) %*% white_y)
    white_resids <- white_y - white_X %*% beta
  }

  nll2 <- n * log(2 * pi) + n + log_det_cormat + n * log(mean(white_resids^2))
  #return(n*log(mean(white_resids^2)))
  return(nll2 / n)
}


#  spatial_cor = (dphi/2)^eta[2]*2*besselK(dphi,eta[2])/gamma(eta[2])

#' This function corresponds to cor_mat function
#' Docs (written by Isaac not Dr. Chu and disclaimer that I might make mistake.
#' Btw, if Dr. Chu wants to use this, please take it.)
#'
#' @param D, this is the distance matrix
#' @param eta, this is the vector of parameters for the covariance function.
#'    Specifically, it should the way should be arranged is,
#'    \itemize{
#'      \item "Exp" - [phi (range parameters), sigma^2 (if use nuggest, then tau^2 / (sigma^2 + tau^2))]
#'      \item ...
#'    }
#' @param cov_model The covariance model, for exponential class (`"Exp"`)
cor.mat <- function(D, eta, cov_model, nug) {
  dphi <- D / eta[1]
  if (cov_model == "Cau") {
    spatial_cor <- (1 + dphi^2)^{
      -eta[2]
    }
  } else if (cov_model == "Exp") {
    spatial_cor <- exp(-dphi)
  } else if (cov_model == "Mat") {
    dtmp <- sqrt(2 * eta[2]) * dphi
    spatial_cor <- 2^(1 - eta[2]) /
      gamma(eta[2]) *
      (dtmp)^eta[2] *
      besselK(dtmp, eta[2])
  } else if (cov_model == "Mat32") {
    spatial_cor <- (1 + 3^.5 * dphi) * exp(-3^0.5 * dphi)
  } else if (cov_model == "Mat52") {
    spatial_cor <- (1 + 5^.5 * dphi + 5 / 3 * dphi^2) * exp(-5^.5 * dphi)
  } else if (cov_model == "Gau") {
    spatial_cor <- exp(-dphi^2)
  } else if (cov_model == "Sph") {
    spatial_cor <- (1 - 1.5 * dphi + 0.5 * dphi^3) * (dphi < 1)
  } else if (cov_model == "GW02") {
    spatial_cor <- (1 - dphi)^2 * (dphi < 1)
  } else {
    print("The input covariance function is not supported")
  }

  diag(spatial_cor) <- 1
  if (nug == FALSE) {
    cormat <- spatial_cor
  } else if (nug == TRUE) {
    q <- length(eta) + 1
    c <- eta[q - 1]
    cormat <- (1 - c) * spatial_cor
  }
  diag(cormat) <- 1 + 1e-8
  return(cormat)
}

# Build Q directly for the "kronecker(Diagonal(m), rep(1/r^2, r^2))" averaging case
build_Q_blockavg <- function(coord, resolution, eta, cov_model, nug = FALSE) {
  stopifnot(nrow(coord) == 2L)
  m <- ncol(coord)
  r <- resolution

  shifts <- -(r - 1):(r - 1) # unique offset differences in 1D
  w1 <- r - abs(shifts) # multiplicities in 1D
  W <- outer(w1, w1, "*") # multiplicities in 2D  => (2r-1)x(2r-1)
  norm <- r**4

  DX <- matrix(rep(shifts, times = length(shifts)), nrow = length(shifts))
  DY <- matrix(rep(shifts, each = length(shifts)), nrow = length(shifts))

  # ---- correlation function (edit to match your BFuns.R / cor.mat definition) ----
  rho <- function(d) {
    cm <- tolower(cov_model)
    if (cm %in% c("exp", "exponential")) {
      # eta[1] = range/scale (example)
      exp(-d / eta[1])
    } else if (cm %in% c("gauss", "gaussian")) {
      exp(-(d / eta[1])^2)
    } else {
      stop("Implement this cov_model in rho(): ", cov_model)
    }
  }

  Q <- matrix(0, m, m)

  for (i in seq_len(m)) {
    Q[i, i] <- 0
    xi <- coord[1, i]
    yi <- coord[2, i]
    for (j in i:m) {
      dx0 <- xi - coord[1, j]
      dy0 <- yi - coord[2, j]

      # Use the W weighting to avoid repeated calculations
      # as the correlation depends solely on the relative dist.
      d <- sqrt((dx0 + DX)^2 + (dy0 + DY)^2) # (2r-1)x(2r-1)
      Qij <- sum(W * rho(d)) / norm

      Q[i, j] <- Qij
      Q[j, i] <- Qij
    }
  }

  # Nugget handling depends on how your cor.mat parameterises it.
  # If your fine-level model is: (1-g)*rho(d) + g*I, then g affects only the diagonal of the fine matrix,
  # and after averaging its diagonal contribution is scaled by 1/r^2.
  #
  # So if eta contains g as the last element when nug=TRUE, a common adjustment is:
  if (isTRUE(nug)) {
    g <- eta[length(eta)] # <-- adjust if your eta layout differs
    diag(Q) <- diag(Q) + g / (r^2)
  }

  Q
}
