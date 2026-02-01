# NOTE: The the numbers here are for MLE estimation. Run the code
#       with the `_reml` suffix data file to get the results for REML.
#       Overall, I can actually see a slightly higher variance (0.438 vs 0.404),
#       but a much lower bias (0.0301 vs 0.0683) for phi. Similar for sigma^2,
#       where bias is now -0.00414 while the variance increased to 0.032.
#       However, both has a slightly higher RMSE than MLE.

# multi_resolution <- readRDS(file = "data/simulation_res.rds")
# multi_resolution <- readRDS(file = "data/simulation_res_reml.rds")
multi_resolution <- readRDS(file = "data/simulation_res_reml_1000.rds")

all(multi_resolution$converge == 1)

# Average likelihood
exp(-mean(multi_resolution$nll / multi_resolution$npoints))
# [1] 0.3271074

nsim <- multi_resolution$nsim

# Average phi, sigma^2
colMeans(multi_resolution$theta[1:nsim, ])
# [1] 4.931721 1.971616
#     phi      sigma^2
mean(multi_resolution$alpha0[1:nsim])
# [1] 2.006571

# RMSE for the estimates.
library(caret)
RMSE(multi_resolution$true_alpha0, multi_resolution$alpha0[1:nsim])
# [1] 0.1671237
RMSE(multi_resolution$true_phi, multi_resolution$theta[1:nsim, 1])
# [1] 0.6383653
RMSE(multi_resolution$true_sigmasq, multi_resolution$theta[1:nsim, 2])
# [1] 0.1761742

sqrt(mean((multi_resolution$alpha0[1:nsim] - multi_resolution$true_alpha0)**2))
# [1] 0.1671237
sqrt(mean((multi_resolution$true_phi - multi_resolution$theta[1:nsim, 1])**2))
# [1] 0.6383653
sqrt(mean((multi_resolution$true_sigmasq - multi_resolution$theta[1:nsim, 2])**2))
# [1] 0.1761742

# Looking Deeper for phi -------------------------------------------------------
hist(x = multi_resolution$theta[1:nsim, 1])

# Bias
mean(multi_resolution$theta[1:nsim, 1] - multi_resolution$true_phi)
# [1] -0.06827891
mean(multi_resolution$theta[1:nsim, 2] - multi_resolution$true_sigmasq)
# [1] -0.02838421
mean(multi_resolution$alpha0[1:nsim] - multi_resolution$true_alpha0)
# [1] 0.006571065

# Variance
var(multi_resolution$theta[1:nsim, 1])
# [1] 0.4036555
var(multi_resolution$theta[1:nsim, 2])
# [1] 0.03029228
var(multi_resolution$alpha0[1:nsim])
# [1] 0.02794304

# Standard Deviation
sqrt(var(multi_resolution$theta[1:nsim, 1]))
# [1] 0.6353389
sqrt(var(multi_resolution$theta[1:nsim, 2]))
# [1] 0.1740468
sqrt(var(multi_resolution$alpha0[1:nsim]))
# [1] 0.1671617

# sqrt(bias**2 + var (or sd**2)) # This should be the RMSE
