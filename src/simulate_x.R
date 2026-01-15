# This file also includes the previous BFuns.R
source(file = "src/likelihood.R", verbose = TRUE)

multi_resolution <- readRDS(file = "data/x_simulation.rds")
slice_idx <- 1

names(multi_resolution)

# The observation values and the corresponding
# coordinate should have the same length
stopifnot(length(multi_resolution$x_obs[slice_idx, ]) == dim(multi_resolution$x.coord)[2])

# Compute the distance matrix.-------------------------------------------------
# Note that given the model, we need to have a distance
# matrix for all the surrounding point.
half <- multi_resolution$resolution %/% 2
offsets <- expand.grid(
  dx = seq(from = -half, to = half, by = 1),
  dy = seq(from = -half, to = half, by = 1)
)
res <- multi_resolution$resolution # Only used for short hand

coarse_coords <- t(multi_resolution$x.coord[slice_idx, , ])
m <- ncol(coarse_coords)

fine_coords <- matrix(NA, nrow = 2, ncol = (res**2) * m)

for (i in seq_len(m)) {
  base <- coarse_coords[, i]
  idx <- ((res**2) * (i - 1) + 1):((res**2) * i)
  fine_coords[, idx] <- base + t(offsets)
}
rm(half, offsets, m, res)
cat(sprintf(
  "The following should form %dx%d block\n",
  multi_resolution$resolution,
  multi_resolution$resolution
))
fine_coords[, 1:9]

D <- as.matrix(dist(t(fine_coords)))
#-------------------------------------------------------------------------------

# Maximimum Likelihood Fit using the -------------------------------------------
lo.bound <- c(0.00001)
up.bound <- c(10)
xx <- MLE.fit(
  multi_resolution$x_obs[slice_idx, ],
  D,
  "Exp",
  c(3),
  nug = FALSE,
  "LB",
  lo.bound,
  up.bound,
  resolution = 3
)
str(xx)
#-------------------------------------------------------------------------------
