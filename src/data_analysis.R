source(file = "src/likelihood.R", verbose = TRUE)
source(file = "src/restricted_likelihood.R", verbose = TRUE)
source(file = "src/kriging.R", verbose = TRUE)

gal <- readRDS(file = "lib_data/ASTData/galData.Rds")
mgs <- readRDS(file = "lib_data/ASTData/mgsData.Rds")

BAU_size <- 0.02
gal_res <- 4 # target resolution
mgs_res <- 9 # src resolution

# Scale the coordinate
gal$X <- gal$X / BAU_size
gal$Y <- gal$Y / BAU_size
mgs$X <- mgs$X / BAU_size
mgs$Y <- mgs$Y / BAU_size

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
# List of 7
#  $ theta             : num [1:2] 10 1043
#  $ converge (1 = yes): num 1
#  $ nll               : num 4330
#  $ alpha0            : num [1, 1] 11.6
#  $ alpha_0_var       : num 2
#  $ eta               : num 10
#  $ Q                 : num [1:1000, 1:1000] 6.43e-01 4.72e-23 5.79e-31 3.37e-35 2.05e-17 ...
# NULL
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
#  0.6428  0.6428  0.6428  0.6428  0.6428  0.6428
eig_min <- min(eigen(Q_hat, symmetric = TRUE, only.values = TRUE)$values)
# [1] 0.0009410232
cond <- kappa(Q_hat)
# [1] 746.0601
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
  sigmasq = xx$theta[2],
  obs.cor = Q_hat
)
saveRDS(kriging_res, file = "data/analysis_kriging_results.rds")
# > str(kriging_res)
# List of 4
#  $ nu      : num [1:869] 0.171 1.1729 0.5383 0.0719 1.995 ...
#  $ lambda  : num [1:869, 1:1000] 7.29e-05 5.00e-04 2.30e-04 3.07e-05 8.51e-04 
# ...
#  $ blup    : num [1:869] 64.5 9.7 11.9 22.2 11.6 ...
#  $ blup_var: num [1:869] 275 747 543 401 859 ...

kriging_res <- readRDS(file = "data/analysis_kriging_results.rds")

# 4. Build Geostatistical Model -----------------------------------------------#

# 5. Some Visualisations for the poster ---------------------------------------#
save_dir <- "../poster/figures/"
pdf(
  paste0(save_dir, "mgs_gal_scatter.pdf"),
  width = 9.0,
  height = 6.0,
  # pointsize = 20
)
par(mar = c(4.5, 4.5, 2.2, 1.0), las = 1)
plot(
  gal$X,
  gal$Y,
  xlim = range(c(mgs$X, gal$X), na.rm = TRUE),
  ylim = range(c(mgs$Y, gal$Y), na.rm = TRUE),
  asp = 1, # preserve geometry
  pch = 16,
  cex = 0.35 * 2,
  col = adjustcolor("red", 0.25),
  xlab = "X (kpc)",
  ylab = "Y (kpc)",
  main = "Locations of Observations of Two Related Astronomic Datasets (MGS & GAL)"
)

grid(col = "grey90")

points(
  mgs$X,
  mgs$Y,
  pch = 16,
  cex = 0.35,
  col = adjustcolor("black", 0.45)
)

legend(
  "topright",
  legend = c("MGS (2x2)", "GAL (82x82)"),
  pch = 16,
  col = c(adjustcolor("black", 0.6), adjustcolor("red", 0.6)),
  pt.cex = 0.8,
  bty = "n"
)

dev.off()
