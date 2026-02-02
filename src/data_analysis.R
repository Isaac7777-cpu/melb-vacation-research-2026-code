source(file = "src/likelihood.R", verbose = TRUE)
source(file = "src/restricted_likelihood.R", verbose = TRUE)

gal <- readRDS(file = "lib_data/ASTData/galData.Rds")
mgs <- readRDS(file = "lib_data/ASTData/mgsData.Rds")

gal_res <- 25 # target resolution
mgs_res <- 11 # src resolution

str(gal)
str(mgs)

# 1. Extract the data of interest ---------------------------------------------#

coord <- rbind(mgs$X, mgs$Y)
obs <- mgs$CO21_mgsd

# 2. Fit the model using REML -------------------------------------------------#

mgs_half <- mgs_res %/% 2L
src_offsets <- expand.grid(
  dx = seq(from = -mgs_half, to = mgs_half, by = 1),
  dy = seq(from = -mgs_half, to = mgs_half, by = 1)
)
m <- ncol(coord)

fine_coords <- matrix(NA, nrow = 2, ncol = (mgs_res**2) * m)

for (i in seq_len(m)) {
  base <- coord[, i]
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

xx <- REML.fit.memeff(
  y = obs,
  coord = coord,
  cov_model = "Exp",
  eta.ini = c(3),
  nug = FALSE,
  opt = "LB",
  lo_bound = c(0.00001),
  up_bound = c(10),
  resolution = mgs_res,
  verbose = TRUE
)

str(xx)
# List of 6
#  $ theta             : num [1:2] 1.23 8.43e+05
#  $ converge (1 = yes): num 1
#  $ nll               : num 4817
#  $ alpha0            : num [1, 1] -20.3
#  $ alpha_0_var       : num 6766
#  $ eta               : num 1.23
saveRDS(xx, file = "data/analysis_fitting_results.rds")

# 2.1 Looking at the results... -----------------------------------------------#
#   Clearly, there is onthing that really stands out, which is the estimated
# sigma^2. I will check the following:
#   1. observed value being large itself
#   2. illed condition Q

c(range = range(obs), mean = mean(obs), var = var(obs))
#      range1      range2        mean         var
#    1.531703  512.768677   15.705747 1745.264561
# This seems to be pretty larged, but just to make sure, we will also look at
# the estimated Q matrix and see if it is ill conditioned.

Q_hat <- build_Q_blockavg(
  coord = coord,
  cov_model = "Exp",
  eta = xx$theta[1],
  resolution = mgs_res,
  nug = FALSE
)

diag_summary <- summary(diag(Q_hat))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0598  0.0598  0.0598  0.0598  0.0598  0.0598
eig_min <- min(eigen(Q_hat, symmetric = TRUE, only.values = TRUE)$values)
# [1] 0.0001271381
cond <- kappa(Q_hat)
# [1] 10134.77
# Okay, so the results is showing descent conditioned matrix with descent
# minimum eigenvalues of the matrix and the kappa is also exceptable. This rules
# out the possibility of it being numerical issues. Given that this maybe
# something inherent to the data / setup, let's just ignore it for now.

# 3. Perform Kriging ----------------------------------------------------------#
