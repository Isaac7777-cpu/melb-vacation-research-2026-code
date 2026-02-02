source(file = "lib/BFuns.R", verbose = TRUE)

#' Calculate the variance
cov.exponential <- function(dist, sigmasq, phi, tausq = NA) {
  cov.mat <- sigmasq * exp(-dist / phi)
  if (!is.na(tausq)) {
    diag(cov.mat) <- diag(cov.mat) + tausq
  }
  cov.mat
}

ordinary.kriging <- function(
  request_coord,
  target_res,
  src_res,
  obs,
  coord,
  alpha0,
  sigmasq,
  phi
) {
  # Compute target coordinates
  target_half <- target_res %/% 2L
  target_offsets <- expand.grid(
    dx = seq(-target_half, target_half),
    dy = seq(-target_half, target_half)
  )
  target_size <- dim(request_coord)[1]
  coord.target.box <- array(
    dim = c(target_size, target_res**2, 2),
    dimnames = list(
      point = NULL,
      cell = NULL,
      coord = c("x", "y")
    )
  )
  for (p in seq_len(target_size)) {
    coord.target.box[p, , ] <- t(request_coord[p, ] + t(target_offsets))
  }

  # Get the distance tensor
  half <- src_res %/% 2
  src_offsets <- expand.grid(
    dx = seq(from = -half, to = half, by = 1),
    dy = seq(from = -half, to = half, by = 1)
  )
  coord.src.box <- array(
    # dim = c(sim_res$npoints, src_res**2, 2),
    dim = c(length(obs), src_res**2, 2),
    dimnames = list(
      point = NULL,
      cell = NULL,
      coord = c("x", "y")
    )
  )
  for (p in seq_along(obs)) {
    coord.src.box[p, , ] <- t(coord[, p] + t(src_offsets))
  }

  # Reduce to the distance

  m <- dim(coord.src.box)[1]
  s <- dim(coord.src.box)[2]
  n_targets <- dim(coord.target.box)[1]
  k <- dim(coord.target.box)[2]

  src_xy <- matrix(coord.src.box, ncol = 2) # Flatten the first two dimensions
  tgt_xy <- matrix(coord.target.box, ncol = 2) # Flatten the first two dimensions

  dx_mat <- outer(src_xy[, 1], tgt_xy[, 1], "-")
  dy_mat <- outer(src_xy[, 2], tgt_xy[, 2], "-")

  dx_raw <- array(dx_mat, dim = c(m, s, n_targets, k))
  dy_raw <- array(dy_mat, dim = c(m, s, n_targets, k))

  dx <- aperm(dx_raw, perm = c(1, 3, 4, 2)) # (m,t,k,s)
  dy <- aperm(dy_raw, perm = c(1, 3, 4, 2)) # (m,t,k,s)

  dist <- sqrt(dx**2 + dy**2) # (m,t,k,s)

  # Obtain the covariance vector
  coord.cov <- cov.exponential(
    dist = dist,
    sigmasq = sigmasq,
    phi = phi,
    tausq = NA
  ) # (m,t,k,s)

  # Average over the neighbour dim (k, s)
  c <- apply(coord.cov, c(1, 2), mean) # (m, t)

  # Retrieved the covariance matrix
  # m <- sim_res$npoints
  m <- length(obs)
  stopifnot(
    "The coordinate vector should be of size 2 x npoints" = m == ncol(coord)
  )
  fine_coords <- matrix(NA, nrow = 2, ncol = (src_res**2) * m)
  for (i in seq_len(m)) {
    base <- coord[, i]
    idx <- ((src_res**2) * (i - 1) + 1):((src_res**2) * i)
    fine_coords[, idx] <- base + t(src_offsets)
  }

  cat(sprintf("The following should form %dx%d block\n", src_res, src_res))
  print(fine_coords[, 1:9])

  pairwise_dist <- as.matrix(dist(t(fine_coords)))
  obs.cov.full <- sigmasq *
    cor.mat(
      D = pairwise_dist,
      cov_model = "Exp",
      eta = c(phi, sigmasq),
      nug = FALSE
    )
  resolution_matrix <- kronecker(
    diag(m),
    matrix(rep(1 / src_res**2, src_res**2), nrow = 1)
  )
  obs.cov <- resolution_matrix %*% obs.cov.full %*% t(resolution_matrix)

  # Compute the BLUP
  L <- t(chol(obs.cov))
  r <- obs - rep(alpha0, m)
  quad_term <- backsolve(r = t(L), x = forwardsolve(l = L, x = r))
  blup <- alpha0 + crossprod(c, quad_term)
  # blup <- drop(blup) # Turn to a single scalar

  # Compute the optimiser
  one <- rep(1, m) # m
  Vinv_one <- backsolve(t(L), forwardsolve(L, one)) # m
  Vinv_c <- backsolve(t(L), forwardsolve(L, c)) # (m,t)
  nu <- drop(1 - crossprod(one, Vinv_c)) / drop(crossprod(one, Vinv_one)) # t
  lambda <- Vinv_c + Vinv_one %*% t(nu) # m
  blup_check <- crossprod(lambda, obs)

  # To check for correctness of the code, we can check that lambda^T %*% obs
  # is indeed the BLUP we calculated using the simplified formula above.
  stopifnot(
    "lambda^top %*% obs != BLUP" = all.equal(
      drop(blup),
      drop(blup_check),
      tolerance = 1e-10
    )
  )

  # Check if they are the same
  Vinv_obs <- backsolve(t(L), forwardsolve(L, obs))
  alpha0_gls <- drop(crossprod(one, Vinv_obs) / crossprod(one, Vinv_one))

  stopifnot(
    "The sum of lambda should be 1 to satisfy the constraint" = all.equal(
      drop(colSums(lambda)),
      rep(1, n_targets),
      tolerance = 1e-10
    )
  )

  stopifnot(
    "Recalculated alpha0 should be the same as the one loaded" = all.equal(
      unname(alpha0),
      alpha0_gls,
      tolerance = 1e-10
    )
  )

  # Compute the BLUP variance
  k <- dim(coord.target.box)[2]
  w <- rep(1 / k, k)

  D_p <- as.matrix(dist(as.matrix(target_offsets)))
  Sigma_p <- cov.exponential(
    dist = D_p,
    phi = phi,
    sigmasq = sigmasq,
    tausq = NA
  )
  block_var <- drop(crossprod(w, Sigma_p %*% w))

  # term2 <- diag(t(c) %*% Vinv_c)
  term2 <- colSums(c * Vinv_c)

  blup_var <- block_var - term2 + drop(crossprod(one, Vinv_one)) * (nu**2)

  list(
    nu = nu,
    lambda = lambda,
    blup = drop(blup),
    blup_var = drop(blup_var)
  )
}
