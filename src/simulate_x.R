source(file = "src/likelihood.R", verbose = TRUE)
source(file = "src/restricted_likelihood.R", verbose = TRUE)

library(optparse)
library(foreach)
library(doParallel)

# Simulation Function ---------------------------------------------------------
simulate_slice <- function(
  use_reml,
  slice_idx,
  res,
  obs,
  coord,
  lo.bound,
  up.bound,
  verbose = FALSE
) {
  # Compute the distance matrix.
  # Note that given the model, we need to have a distance
  # matrix for all the surrounding point.
  half <- res %/% 2
  offsets <- expand.grid(
    dx = seq(from = -half, to = half, by = 1),
    dy = seq(from = -half, to = half, by = 1)
  )

  m <- ncol(coord)

  fine_coords <- matrix(NA, nrow = 2, ncol = (res**2) * m)

  for (i in seq_len(m)) {
    base <- coord[, i]
    idx <- ((res**2) * (i - 1) + 1):((res**2) * i)
    fine_coords[, idx] <- base + t(offsets)
  }

  if (verbose) {
    cat(sprintf("The following should form %dx%d block\n", res, res))
    print(fine_coords[, 1:9])
  }

  D <- as.matrix(dist(t(fine_coords)))

  # Maximimum Likelihood Fit using the
  if (use_reml) {
    xx <- REML.fit(
      obs,
      D,
      "Exp",
      c(3),
      nug = FALSE,
      "LB",
      lo.bound,
      up.bound,
      resolution = res,
      verbose = verbose
    )
  } else {
    xx <- MLE.fit(
      obs,
      D,
      "Exp",
      c(3),
      nug = FALSE,
      "LB",
      lo.bound,
      up.bound,
      resolution = res,
      verbose = verbose,
      nll_correction = profile.nll.correction
    )
  }

  if (verbose) {
    cat(str(xx))
  }

  xx
}

# Scripts run for everything exists in the loaded dataset.

option_list <- list(
  make_option(
    c("-v", "--verbose"),
    type = "logical",
    default = FALSE,
    help = "Print extra output [default= %default]"
  ),
  make_option(
    c("-n", "--num-sim"),
    type = "integer",
    default = 100,
    help = "Number of simulations to run, [default= %default (note that this can depends on the data.)]",
    metavar = "INTEGER"
  ),
  make_option(
    c("-c", "--num-core"),
    type = "integer",
    default = 2,
    help = "Number of core to run the simulations.",
    metavar = "INTEGER"
  ),
  make_option(
    c("-d", "--data-file"),
    type = "character",
    default = "data/x_simulation.rds",
    help = "The simulation data file [default= %default]"
  ),
  make_option(
    c("-r", "--use-reml"),
    type = "logical",
    default = FALSE,
    help = "Whether to use REML instead of MLE"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

multi_resolution <- readRDS(file = opt[["data-file"]])
verbose <- opt[["verbose"]]
num_sim <- opt[["num-sim"]]
use_reml <- opt[["use-reml"]]

stopifnot(num_sim <= dim(multi_resolution$x_obs)[1])

if (interactive()) {
  str(multi_resolution)
  num_sim <- 2
  verbose <- FALSE
}

sim_res <- list(
  # Simulation settings
  seed = multi_resolution$seed,
  resolution = multi_resolution$resolution,
  npoints = multi_resolution$npoints,
  nsim = num_sim,
  nug = multi_resolution$nug,
  use_reml = use_reml,
  # True parameters
  true_tausq = multi_resolution$tausq,
  true_alpha0 = multi_resolution$alpha0,
  true_sigmasq = multi_resolution$sigmasq,
  true_phi = multi_resolution$phi,
  # Estimated Parameters
  alpha0 = rep(0.0, num_sim),
  alpha0_var = rep(0.0, num_sim),
  theta = matrix(data = NA, nrow = num_sim, ncol = 2),
  converge = rep(0.0, num_sim),
  nll = rep(0.0, num_sim)
)

set.seed(multi_resolution$seed)

cl <- makeCluster(opt[["num-core"]])
registerDoParallel(opt$"num-core")

res_list <- foreach(
  slice_idx = seq_len(num_sim),
  .packages = character(0)
) %dopar%
  {
    xx <- simulate_slice(
      use_reml = use_reml,
      slice_idx = slice_idx,
      res = multi_resolution$resolution,
      obs = multi_resolution$x_obs[slice_idx, ],
      coord = multi_resolution$x.coord[slice_idx, , ],
      lo.bound = c(0.00001),
      up.bound = c(10),
      verbose = verbose
    )

    list(
      alpha0 = xx$alpha0,
      alpha0_var = xx$alpha_0_var,
      theta = xx$theta,
      converge = xx[["converge (1 = yes)"]],
      nll = xx$nll
    )
  }

stopCluster(cl)

# Combine
sim_res$alpha0 <- vapply(res_list, `[[`, numeric(1), "alpha0")
sim_res$alpha0_var <- vapply(res_list, `[[`, numeric(1), "alpha0_var")
sim_res$converge <- vapply(res_list, `[[`, numeric(1), "converge")
sim_res$nll <- vapply(res_list, `[[`, numeric(1), "nll")
sim_res$theta <- do.call(rbind, lapply(res_list, `[[`, "theta"))

if (use_reml) {
  saveRDS(object = sim_res, file = "data/simulation_res_reml.rds")
} else {
  saveRDS(object = sim_res, file = "data/simulation_res.rds")
}

if (interactive()) {
  load_res <- readRDS("data/simulation_res.rds")
  load_res_1 <- readRDS("data/simulation_res.rds")

  colMeans(x = load_res$theta)
  mean(x = load_res$alpha0)
}
