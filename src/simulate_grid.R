library(fields)
source("lib/BFuns.R")

set.seed(47)

n <- 100
ndraw <- 500
grid <- list(
  x = seq(from = 0, by = 1, length.out = n),
  y = seq(from = 0, by = 1, length.out = n)
)
coords <- as.matrix(expand.grid(grid))

dist <- as.matrix(dist(coords))

nug <- FALSE
tausq <- 0.5
sigmasq <- 4 / 2
phi <- 5

# If we consider nugget effect, we need to correct the variance used to
# generate the correlation be corrected as
#   (tau^2) / (sigma^2 + tau^2)
# This result I (Isaac Leong) believe comes from simply algebraic simplification.
if (nug) {
  var_nug_corrected <- tausq / (sigmasq + tausq)
} else {
  var_nug_corrected <- sigmasq
}

# Call the cor.mat from `"BFuns.R"`.
# *** BE CAREFUL OF THE COVARIANCE PARAMETER LAYOUT ***
cormat <- cor.mat(
  D = dist,
  eta = c(phi, var_nug_corrected),
  cov_model = "Exp",
  nug = nug
)

# If using nuggest, we use the equation:
#    var-covariance-mat = (sigma^2 + tau^2) x R (=sigma^2 x R + tau^2 x I)
# to obtain the variance covariance matrix.
# If we are not considering nugget effect, we then do,
#    var-covariance-mat = sigma^2 x R
if (nug) {
  var_covar <- (sigmasq + tausq) * cormat
} else {
  var_covar <- sigmasq * cormat
}

#  Generate the grid but flatten
ymat_flat <- rmvnorm(n = ndraw, sigma = var_covar)

ymat <- aperm(
  array(as.vector(t(ymat_flat)), dim = c(100, 100, ndraw)),
  c(3, 1, 2)
)

if (nug) {
  save(ymat, file = sprintf("data/multi-res-%ddraws.RData", ndraw))
} else {
  save(ymat, file = sprintf("data/multi-res-%ddraws_no_nug.RData", ndraw))
}

# load(file = "ymat_10000draws.RData")
