get_sample_from_slice <- function(
  ymat,
  slice_idx,
  seed,
  alpha0,
  npoints,
  resolution,
  verbose = FALSE
) {
  set.seed(seed = seed)

  # Obtain the basic simulation field for X(v)
  eta_field <- ymat[slice_idx, , ]

  # In order to get the actual field of X(v), we need add an expectation
  x_field <- eta_field + alpha0

  # Picking points in the random field.
  sample_width <- dim(x_field)[1] - resolution + 1
  sample_height <- dim(x_field)[2] - resolution + 1
  sample_length <- sample_width * sample_height

  # Sample the index for the "picked" x points
  x.index <- sample.int(n = sample_length, size = npoints, replace = FALSE)

  # Reorganised as list of coordinates
  x.coord <- sapply(x.index, function(ind) {
    c(
      ((ind - 1) %% sample_width) + ((resolution + 1) %/% 2),
      ((ind - 1) %/% sample_height) + ((resolution + 1) %/% 2)
    )
  })

  if (verbose) {
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
  }

  # Now, after checking that all the boxes are valid,
  # we can get the samples at the desired resolution.
  half <- resolution %/% 2
  x_obs <- sapply(seq_len(npoints), function(k) {
    cx <- x.coord[1, k]
    cy <- x.coord[2, k]
    mean(x_field[(cx - half):(cx + half), (cy - half):(cy + half)])
  })

  if (verbose) {
    length(x_obs)
  }

  return(list(
    x_obs = x_obs,
    x.coord = x.coord
  ))
}

# Interactive Testing Code Chunk------------------------------------------------
if (interactive()) {
  # Pick whether to use nugget effect for the  base random field to use.
  load(file = "data/multi-res-500draws_no_nug.RData", verbose = TRUE)
  # load(file = "data/multi-res-10000draws.RData", verbose = TRUE)
  # This loads 10000 independent draw for an exponential
  # covariance gaussian random field of size 100x100.
  # The variable is called ymat.

  # Simulation a slide to use
  seed <- 10
  alpha0 <- 2.0
  npoints <- 400
  nslice <- dim(ymat)[1]
  resolution <- 3

  sample_obj <- list(
    seed = seed,
    alpha0 = alpha0,
    npoints = npoints,
    resolution = resolution,
    x_obs = array(dim = c(nslice, npoints)),
    x.coord = array(dim = c(nslice, 2, npoints))
  )

  for (slice_idx in seq_len(nslice)) {
    sample <- get_sample_from_slice(
      ymat = ymat,
      slice_idx = slice_idx,
      seed = seed,
      alpha0 = alpha0,
      npoints = npoints,
      resolution = resolution,
      verbose = TRUE
    )

    sample_obj[["x_obs"]][slice_idx, ] <- sample$x_obs
    sample_obj[["x.coord"]][slice_idx, , ] <- sample$x.coord
  }

  saveRDS(
    sample_obj,
    file = "data/x_simulation.rds"
  )

  obs_object <- readRDS(file = "data/x_bundle.rds")
}
