# Pick whether to use nugget effect for the  base random field to use.
load(file = "data/multi-res-100draws_no_nug.RData", verbose = TRUE)
# load(file = "data/multi-res-10000draws.RData", verbose = TRUE)
# This loads 10000 independent draw for an exponential
# covariance gaussian random field of size 100x100.
# The variable is called ymat.

seed <- 1
set.seed(seed)

# Simulation slide to use
slice_idx <- 1

# Obtain the basic simulation field for X(v)
eta_field <- ymat[1, , ]

# In order to get the actual field of X(v), we need add an expectation
alpha0 <- 2.0
x_field <- eta_field + alpha0

# Picking points in the random field.
ndraws <- 400
resolution <- 3
sample_width <- dim(x_field)[1] - resolution + 1
sample_height <- dim(x_field)[2] - resolution + 1
sample_length <- sample_width * sample_height

# Sample the index for the "picked" x points
x.index <- sample.int(n = sample_length, size = ndraws, replace = FALSE)

# Reorganised as list of coordinates
x.coord <- sapply(x.index, function(ind) {
  c(
    (ind %/% sample_width) + (resolution %/% 2 + 1),
    (ind %% sample_height) + (resolution %/% 2 + 1)
  )
})

# Just a quick check of the generated points
# Everything should be integer and between the range of 100 x 100
plot(
  x = x.coord[1, ],
  y = x.coord[2, ],
  main = "Should be bounded in the Domain",
  pch = 16,
  cex = 0.3,
  xlab = "x",
  ylab = "y"
)
# Drawing 3x3 box
rect(
  xleft = x.coord[1, ] - 1,
  xright = x.coord[1, ] + 1,
  ybottom = x.coord[2, ] - 1,
  ytop = x.coord[2, ] + 1
)
abline(v = c(1, 100), col = "red")
abline(h = c(1, 100), col = "red")

# Now, after checking that all the boxes are valid,
# we can get the samples at the desired resolution.
half <- resolution %/% 2
x_obs <- sapply(seq_len(ndraws), function(k) {
  cx <- x.coord[1, k]
  cy <- x.coord[2, k]
  mean(x_field[(cx - half):(cx + half), (cy - half):(cy + half)])
})
length(x_obs)

saveRDS(
  list(
    resolution = resolution,
    x_obs = x_obs,
    x.coord = x.coord,
    slice_idx = slice_idx,
    alpha0 = alpha0,
    seed = seed
  ),
  file = "data/x_bundle.rds"
)

obs_object <- readRDS(file = "data/x_bundle.rds")
