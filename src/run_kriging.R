# WARNING: This script assume the simulation result and the
#          data has the same slice index
source(file = "src/kriging.R", verbose = TRUE)
sim_res <- readRDS(file = "data/simulation_res_reml.rds")
# multi_resolution <- readRDS(file = "data/x_simulation.rds")

slice_idx <- 2
target_res <- 5
request_coord <- c(x = 14, y = 92)
# If you want random coordinate for the 100x100 grid domain
# request_coord <- sample.int(size = 2, n = 98, replace = TRUE) + 1

alpha0 <- sim_res$alpha0[slice_idx]
phi <- sim_res$theta[slice_idx, 1]
sigmasq <- sim_res$theta[slice_idx, 2]
coord <- sim_res$x.coord[slice_idx, , ]
obs <- sim_res$x_obs[slice_idx, ]

# Get the distance tensor
src_res <- sim_res$resolution
half <- src_res %/% 2
src_offsets <- expand.grid(
  dx = seq(from = -half, to = half, by = 1),
  dy = seq(from = -half, to = half, by = 1)
)
coord.box <- array(
  dim = c(sim_res$npoints, src_res**2, 2),
  dimnames = list(
    point = NULL,
    cell = NULL,
    coord = c("x", "y")
  )
)
for (p in seq_len(400)) {
  coord.box[p, , ] <- t(coord[, p] + t(src_offsets))
}

# Reduce to the distance
dx <- coord.box[,, "x", drop = FALSE] - request_coord[1]
dy <- coord.box[,, "y", drop = FALSE] - request_coord[2]
coord.dist <- sqrt(dx**2 + dy**2)

# Obtain the covariance vector
coord.cov <- cov.exponential(
  dist = coord.dist,
  sigmasq = sigmasq,
  phi = phi,
  tausq = NA
)
c <- rowMeans(x = coord.cov)

# Retrieved the covariance matrix
m <- sim_res$npoints
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
obs.cov.full <- cor.mat(
  D = pairwise_dist,
  cov_model = "Exp",
  eta = c(phi),
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
one <- rep(1, m)
Vinv_c <- backsolve(t(L), forwardsolve(L, c))
Vinv_one <- backsolve(t(L), forwardsolve(L, one))
nu <- drop((1 - crossprod(one, Vinv_c)) / crossprod(one, Vinv_one))
lambda <- Vinv_c + nu * Vinv_one
# Check if they are the same
Vinv_obs <- backsolve(t(L), forwardsolve(L, obs))
alpha0_gls <- drop(crossprod(one, Vinv_obs) / crossprod(one, Vinv_one))

stopifnot(
  "The sum of lambda should be 1 to satisfy the constraint" = all.equal(
    sum(lambda),
    1.0,
    tolerance = 1e-10
  )
)

alpha0 - alpha0_gls

stopifnot(
  "Recalculated alpha0 should be the same as the one loaded" = all.equal(
    unname(alpha0),
    alpha0_gls,
    tolerance = 1e-10
  )
)

# Compute the BLUP variance
blup_var <- drop(
  cov.exponential(
    dist = 0,
    phi = phi,
    sigmasq = sigmasq,
    tausq = NA
  ) -
    crossprod(c, Vinv_c) +
    nu**2 * crossprod(one, Vinv_one)
)
