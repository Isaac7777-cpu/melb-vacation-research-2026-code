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
