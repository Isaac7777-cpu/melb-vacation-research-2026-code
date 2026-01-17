multi_resolution <- readRDS(file = "data/simulation_res.rds")

all(multi_resolution$converge == 1)

# Average alpha0
mean(multi_resolution$alpha0)

# Average phi, sigma^2
colMeans(multi_resolution$theta[1:500, ])
# [1] 4.931721 1.971616
#     phi      sigma^2

# RMSE for the estimates.
library(caret)
RMSE(multi_resolution$true_alpha0, multi_resolution$alpha0)
# [1] 0.1671237
RMSE(multi_resolution$true_phi, multi_resolution$theta[, 1])
# [1] 0.6383653
RMSE(multi_resolution$true_sigmasq, multi_resolution$theta[, 2])
# [1] 0.1761742

sqrt(mean((multi_resolution$alpha0 - multi_resolution$true_alpha0)**2))
# [1] 0.1671237
sqrt(mean((multi_resolution$true_phi - multi_resolution$theta[, 1])**2))
# [1] 0.6383653
sqrt(mean((multi_resolution$true_sigmasq - multi_resolution$theta[, 2])**2))
# [1] 0.1761742

# Looking Deeper for phi -------------------------------------------------------
hist(x = multi_resolution$theta[, 1])

# Bias
bias <- mean(multi_resolution$theta[, 1] - multi_resolution$true_phi)
# [1] -0.06827891
mean(multi_resolution$theta[, 2] - multi_resolution$true_sigmasq)
# [1] -0.02838421
mean(multi_resolution$alpha0 - multi_resolution$true_alpha0)
# [1] 0.006571065

# Variance
var <- var(multi_resolution$theta[, 1])
# [1] 0.4036555
var(multi_resolution$theta[, 2])
# [1] 0.03029228
var(multi_resolution$alpha0)
# [1] 0.02794304

# sqrt(bias**2 + var) # This should be the MSE
