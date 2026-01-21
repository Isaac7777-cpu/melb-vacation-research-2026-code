# WARNING: This script assume the simulation result and the
#          data has the same slice index
source(file = "src/kriging.R", verbose = TRUE)
sim_res <- readRDS(file = "data/simulation_res.rds")
multi_resolution <- readRDS(file = "data/x_simulation.rds")

slice_idx <- 1
request_coord <- c(x = 14, y = 92)
# If you want random coordinate for the 100x100 grid domain
# request_coord <- sample.int(size = 2, n = 98, replace = TRUE) + 1

alpha0 <- sim_res$alpha0[slice_idx]
phi <- sim_res$theta[slice_idx, 1]
sigmasq <- sim_res$theta[slice_idx, 2]
coord <- multi_resolution$x.coord[slice_idx, , ]
obs <- multi_resolution$x_obs[slice_idx, ]

# Get the distance tensor
res <- multi_resolution$resolution
half <- res %/% 2
offsets <- expand.grid(
  dx = seq(from = -half, to = half, by = 1),
  dy = seq(from = -half, to = half, by = 1)
)
coord.box <- array(
  dim = c(multi_resolution$npoints, res**2, 2),
  dimnames = list(
    point = NULL,
    cell = NULL,
    coord = c("x", "y")
  )
)
for (p in seq_len(400)) {
  coord.box[p, , ] <- t(coord[, p] + t(offsets))
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

