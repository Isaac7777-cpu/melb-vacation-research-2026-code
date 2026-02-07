# mc_table_plots.R

library(tidyverse)
library(scales)

mc <- tribble(
  ~Param    , ~Method , ~R   , ~Mean  , ~Bias   , ~SD    , ~RMSE  ,

  # Table 1
  "sigma^2" , "MLE"   ,  100 , 1.9710 , -0.0290 , 0.1975 , 0.1987 ,
  "sigma^2" , "REML"  ,  100 , 1.9955 , -0.0045 , 0.2040 , 0.2030 ,
  "sigma^2" , "MLE"   ,  300 , 1.9644 , -0.0356 , 0.1851 , 0.1882 ,
  "sigma^2" , "REML"  ,  300 , 1.9888 ,  0.0379 , 0.1719 , 0.1912 ,
  "sigma^2" , "MLE"   ,  500 , 1.9716 , -0.0284 , 0.1740 , 0.1762 ,
  "sigma^2" , "REML"  ,  500 , 1.9959 , -0.0041 , 0.1797 , 0.1796 ,

  "phi"     , "MLE"   ,  100 , 4.9624 , -0.0376 , 0.6541 , 0.6519 ,
  "phi"     , "REML"  ,  100 , 5.0617 ,  0.0617 , 0.6810 , 0.6804 ,
  "phi"     , "MLE"   ,  300 , 4.9388 , -0.0612 , 0.6576 , 0.6594 ,
  "phi"     , "REML"  ,  300 , 5.0379 , -0.0112 , 0.6854 , 0.6853 ,
  "phi"     , "MLE"   ,  500 , 4.9317 , -0.0683 , 0.6353 , 0.6384 ,
  "phi"     , "REML"  ,  500 , 5.0301 ,  0.0301 , 0.6623 , 0.6622 ,

  "alpha_0" , "MLE"   ,  100 , 1.9921 , -0.0079 , 0.1558 , 0.1552 ,
  "alpha_0" , "REML"  ,  100 , 1.9922 , -0.0078 , 0.1558 , 0.1552 ,
  "alpha_0" , "MLE"   ,  300 , 2.0074 ,  0.0074 , 0.1719 , 0.1718 ,
  "alpha_0" , "REML"  ,  300 , 2.0075 ,  0.0075 , 0.1719 , 0.1718 ,
  "alpha_0" , "MLE"   ,  500 , 2.0066 ,  0.0066 , 0.1672 , 0.1671 ,
  "alpha_0" , "REML"  ,  500 , 2.0065 ,  0.0065 , 0.1672 , 0.1671 ,

  # Table 2 (REML only)
  "sigma^2" , "REML"  ,  600 , 1.9950 , -0.0050 , 0.1968 , 0.1967 ,
  "sigma^2" , "REML"  ,  800 , 1.9990 , -0.0188 , 0.1978 , 0.1977 ,
  "sigma^2" , "REML"  , 1000 , 1.9978 , -0.0022 , 0.1977 , 0.1976 ,

  "phi"     , "REML"  ,  600 , 4.9999 , -0.0001 , 0.6878 , 0.6872 ,
  "phi"     , "REML"  ,  800 , 5.0189 , -0.0010 , 0.7029 , 0.7027 ,
  "phi"     , "REML"  , 1000 , 5.0165 ,  0.0165 , 0.6976 , 0.6974 ,

  "alpha_0" , "REML"  ,  600 , 2.0122 ,  0.0122 , 0.1784 , 0.1967 ,
  "alpha_0" , "REML"  ,  800 , 2.0065 ,  0.0065 , 0.1796 , 0.1796 ,
  "alpha_0" , "REML"  , 1000 , 2.0104 ,  0.0104 , 0.1784 , 0.1786
)

## Base-R poster plots (start from your `mc <- tribble(...)`)
## Uses default x-axis ticks/labels (so no forced tick at 600).

param_levels <- c("sigma^2", "phi", "alpha_0")
param_labels <- c(expression(sigma^2), expression(phi), expression(alpha[0]))
method_levels <- c("MLE", "REML")

mc$Param <- factor(mc$Param, levels = param_levels)
mc$Method <- factor(mc$Method, levels = method_levels)

cols <- setNames(c(1, 2, 4), param_levels) # black/red/blue
pch_param <- setNames(c(16, 17, 15), param_levels) # circle/triangle/square

plot_mc_method <- function(
  dat,
  method,
  R_use,
  file = NULL,
  width_in = 12,
  height_in = 4.2,
  point_cex = 1.3,
  lwd = 2
) {
  d <- dat[dat$Method == method & dat$R %in% R_use, , drop = FALSE]
  d <- d[order(d$R), ]

  if (!is.null(file)) {
    pdf(file, width = width_in, height = height_in, pointsize = 18)
    on.exit(dev.off(), add = TRUE)
  }

  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)

  # par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1), oma = c(0, 0, 2, 0))
  par(
    mfrow = c(1, 3),
    mar = c(3.0, 3.0, 1.1, 0.6),
    oma = c(0, 0, 2.0, 0),
    mgp = c(1.8, 0.7, 0)
  )

  metrics <- c("Bias", "SD", "RMSE")

  for (met in metrics) {
    yr <- range(d[[met]], na.rm = TRUE)
    pad <- 0.08 * diff(yr)
    if (!is.finite(pad) || pad == 0) {
      pad <- 0.05
    }
    ylim <- c(yr[1] - pad, yr[2] + pad)

    plot(NA, xlim = range(d$R), ylim = ylim, xlab = "R", ylab = met) # <-- default x-axis ticks/labels

    if (met == "Bias") {
      abline(h = 0, lty = 2)
    }

    for (p in param_levels) {
      dp <- d[d$Param == p, , drop = FALSE]
      lines(dp$R, dp[[met]], col = cols[p], lwd = lwd)
      points(
        dp$R,
        dp[[met]],
        col = cols[p],
        pch = pch_param[p],
        cex = point_cex
      )
    }
  }

  mtext(
    paste0(method, " Monte Carlo summaries"),
    outer = TRUE,
    cex = 1.8,
    font = 2
  )

  par(xpd = NA)
  legend(
    "topright",
    legend = param_labels,
    col = cols[param_levels],
    pch = pch_param[param_levels],
    lwd = lwd,
    bty = "n",
    horiz = TRUE,
    inset = c(0, -0.15)
  )
  par(xpd = FALSE)

  invisible(NULL)
}

# Export poster-ready PDFs
plot_mc_method(
  mc,
  "MLE",
  R_use = c(100, 300, 500),
  file = "fig_mle_base.pdf",
  width_in = 12,
  height_in = 4.2
)

plot_mc_method(
  mc,
  "REML",
  R_use = c(100, 300, 500, 600, 800, 1000),
  file = "fig_reml_base.pdf",
  width_in = 12,
  height_in = 4.2
)

# If you want even bigger for the poster, just bump width/height, e.g.:
# plot_mc_method(mc, "REML", c(100,300,500,600,800,1000), "fig_reml_base_big.pdf", 16, 5)
