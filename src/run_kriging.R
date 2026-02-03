# WARNING: This script assume the simulation result and the
#          data has the same slice index
source(file = "src/kriging.R", verbose = TRUE)
sim_res <- readRDS(file = "data/simulation_res_reml.rds")
# multi_resolution <- readRDS(file = "data/x_simulation.rds")

slice_idx <- 1
target_res <- as.integer(5)
request_coord <- matrix(data = c(x = c(14, 21), y = c(94, 10)), ncol = 2)

stopifnot(target_res >= 1L, target_res %% 2L == 1L) # odd block size
alpha0 <- sim_res$alpha0[slice_idx]
phi <- sim_res$theta[slice_idx, 1]
sigmasq <- sim_res$theta[slice_idx, 2]
coord <- sim_res$x.coord[slice_idx, , ]
obs <- sim_res$x_obs[slice_idx, ]
src_res <- sim_res$resolution

kriging_res <- ordinary.kriging(
  request_coord = request_coord,
  target_res = target_res,
  src_res = src_res,
  obs = obs,
  coord = coord,
  alpha0 = alpha0,
  sigamsq = sigmasq,
  phi = phi
)

str(kriging_res)
