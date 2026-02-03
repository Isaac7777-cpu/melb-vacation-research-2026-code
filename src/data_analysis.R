source(file = "src/likelihood.R", verbose = TRUE)
source(file = "src/restricted_likelihood.R", verbose = TRUE)
source(file = "src/kriging.R", verbose = TRUE)

gal <- readRDS(file = "lib_data/ASTData/galData.Rds")
mgs <- readRDS(file = "lib_data/ASTData/mgsData.Rds")

gal_res <- 82 # target resolution
mgs_res <- 2 # src resolution

str(gal)
str(mgs)

# 1. Extract the data of interest ---------------------------------------------#

mgs_coord <- rbind(mgs$X, mgs$Y)
mgs_obs <- mgs$CO21_mgsd

# 2. Fit the model using REML -------------------------------------------------#

mgs_half <- mgs_res %/% 2L
parity_correction <- ((mgs_res + 1) %% 2) / 2
src_offsets <- expand.grid(
  dx = seq(
    from = -mgs_half + parity_correction,
    to = mgs_half - parity_correction,
    by = 1
  ),
  dy = seq(
    from = -mgs_half + parity_correction,
    to = mgs_half - parity_correction,
    by = 1
  )
)
m <- ncol(mgs_coord)

fine_coords <- matrix(NA, nrow = 2, ncol = (mgs_res**2) * m)

for (i in seq_len(m)) {
  base <- mgs_coord[, i]
  idx <- ((mgs_res**2) * (i - 1) + 1):((mgs_res**2) * i)
  fine_coords[, idx] <- base + t(src_offsets)
}

if (interactive()) {
  cat(sprintf(
    "The following should form %dx%d block for mgs data\n",
    mgs_res,
    mgs_res
  ))
  print(fine_coords[, 1:mgs_res**2])
}

if (mgs_res < 5) {
  D <- as.matrix(dist(t(fine_coords)))

  xx <- REML.fit(
    y = mgs_obs,
    D = D,
    cov_model = "Exp",
    eta.ini = c(3),
    nug = FALSE,
    "LB",
    lo_bound = c(0.00001),
    up_bound = c(10),
    resolution = mgs_res,
    verbose = TRUE
  )
} else {
  # WARNING: The following only works for odd resolution
  xx <- REML.fit.memeff(
    y = mgs_obs,
    coord = mgs_coord,
    cov_model = "Exp",
    eta.ini = c(3),
    nug = FALSE,
    opt = "LB",
    lo_bound = c(0.00001),
    up_bound = c(10),
    resolution = mgs_res,
    verbose = TRUE
  )
}

str(xx)
# List of 6
#  $ theta             : num [1:2] 0.936 5715.562
#  $ converge (1 = yes): num 1
#  $ nll               : num 4430
#  $ alpha0            : num [1, 1] 2.34
#  $ alpha_0_var       : num 64.2
#  $ eta               : num 0.936
saveRDS(xx, file = "data/analysis_fitting_results.rds")

# 2.1 Looking at the results... -----------------------------------------------#
#   Clearly, there is onthing that really stands out, which is the estimated
# sigma^2. I will check the following:
#   1. observed value being large itself
#   2. illed condition Q

c(range = range(mgs_obs), mean = mean(mgs_obs), var = var(mgs_obs))
#      range1      range2        mean         var
#    1.531703  512.768677   15.705747 1745.264561
# This seems to be pretty larged, but just to make sure, we will also look at
# the estimated Q matrix and see if it is ill conditioned.

Q_hat <- xx$Q

diag_summary <- summary(diag(Q_hat))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.4769  0.4769  0.4769  0.4769  0.4769  0.4769
eig_min <- min(eigen(Q_hat, symmetric = TRUE, only.values = TRUE)$values)
# [1] 0.005282388
cond <- kappa(Q_hat)
# [1] 562.4237
# Okay, so the results is showing descent conditioned matrix with descent
# minimum eigenvalues of the matrix and the kappa is also exceptable. This rules
# out the possibility of it being numerical issues. Given that this maybe
# something inherent to the data / setup, let's just ignore it for now.

# 3. Perform Kriging ----------------------------------------------------------#

y_coord_raw <- cbind(gal$X, gal$Y)
y_obs_raw <- gal$Z_RS32
table(is.nan(y_obs_raw))
table(is.infinite(y_obs_raw))

# 3.1. Missing Y Observation --------------------------------------------------#
faulty <- is.na(y_obs_raw) | is.nan(y_obs_raw) | is.infinite(y_obs_raw)
table(faulty) # Check again

y_obs <- y_obs_raw[!faulty]
y_coord <- y_coord_raw[!faulty, ]

kriging_res <- ordinary.kriging(
  request_coord = y_coord,
  target_res = gal_res,
  src_res = mgs_res,
  obs = mgs_obs,
  coord = mgs_coord,
  alpha0 = drop(xx$alpha0),
  phi = xx$theta[1],
  sigmasq = xx$theta[2]
)
saveRDS(kriging_res, file = "data/analysis_kriging_results.rds")
# The results, nu is quite close to being identical...
# > str(kriging_res)
# List of 4
#  $ nu      : num [1:869] 59.5 59.5 59.5 59.5 59.5 ...
#  $ lambda  : num [1:869, 1:1000] 0.000301 0.00027 0.000186 0.000146 0.000173 ...
#  $ blup    : num [1:869] 2.34 2.36 2.32 2.3 2.35 ...
#  $ blup_var: num [1:869] 59.6 59.6 59.6 59.6 59.6 ...
# Just for interest, this is the result for if we have set the resolution to 5
# > str(test)
# List of 4
#  $ nu      : num [1:869] 8.18 11.45 12.17 6.38 43.67 ...
#  $ lambda  : num [1:869, 1:1000] 2.03e-05 4.24e-05 4.25e-05 2.48e-05 2.14e-04 ...
#  $ blup    : num [1:869] 4.41 15.39 2.74 2.38 2.17 ...
#  $ blup_var: num [1:869] 140.2 178.2 187 97.5 745.3 ...

kriging_res <- readRDS(file = "data/analysis_kriging_results.rds")

# 4. Build Geostatistical Model -----------------------------------------------#
