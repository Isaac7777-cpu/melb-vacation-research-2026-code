# Simulating data in chapter 4 setting

# load(
#   "~/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/A1_Spatial/1-RGeo/locs/loc800.Rdata"
# )
#
# setwd(
#   "~/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/A1_Spatial/1-RGeo"
# )
library(mvtnorm)
source("BFuns.R")

# Hyperparameters
sigmasq <- 4 / 2
tausq <- 0.5
r <- 5
theta0 <- c(r, sigmasq, tausq)
theta.ini <- c(4.4, 0.1)
nug <- TRUE
cov_model <- "Exp"
Srep <- 100

# Creating the distance
n <- loc$n
D <- matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
  D[i, ] <- sqrt((loc$x[i] - loc$x)^2 + (loc$y[i] - loc$y)^2)
}

# Define the covariance matrix
q <- length(theta0)
eta0 <- theta0[-q]
if (nug == TRUE) {
  eta0[q - 1] = theta0[q] / (theta0[q - 1] + theta0[q])
}
cormat = cor.mat(D, eta0, cov_model, nug)

if (nug == TRUE) {
  Sigma0 <- (theta0[q - 1] + theta0[q]) * cormat
} else {
  Sigma0 <- theta0[q] * cormat
}

Ymat <- rmvnorm(Srep, sigma = Sigma0)

# We want 400 observations for the 3x3 resolutions

dd <- (loc$x^2 + loc$y^2)^.5
Zdata <- 2 + 3 * dd + Ymat[1, ]


#library(geoR)
#datafm = data.frame(x = loc$x, y = loc$y, d = dd, Z = Zdata)
#data_geo = as.geodata(datafm, data.col = 4, covar.col = 3, covar.names = "dd")

#plot(variog(data_geo, trend = ~dd, option="bin", max.dist=15),
#     xlab = "h", ylab = "variogram")

#lines.variomodel(seq(0, 900, l = 15),
#                 cov.pars = c(2, 5),
#                 cov.model = "mat", kap = 0.5, nug = 0.5)

#zz = likfit(data_geo, trend = ~dd, ini = c(2,2), nug = 0.3)
#summary(zz)

lo.bound = c(0.00001, 0.00001)
up.bound = c(10, 0.5)
xx = MLE.fit(
  Zdata,
  cbind(1, dd),
  D,
  "Exp",
  c(3, 0.3),
  nug = TRUE,
  "LB",
  lo.bound,
  up.bound
)
