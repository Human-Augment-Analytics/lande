# The three ways adaptive_landscape() averages fitness over a population:
# "simulate" must give what it always gave, "lattice" must follow it draw for
# draw, and "convolve" must agree with the exact answer.

# Bumpus's sparrows as the app fits them, within the sexes
bumpus_surface <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) {
      traits <- c("weight", "total_length")
      d <- bumpus[complete.cases(bumpus[, c("survival", traits)]), ]
      prep <- quiet(prepare_selection_data(d, "survival", traits, standardize = TRUE, group = "sex",
                                           add_relative = TRUE, na_action = "none"))
      surf <- quiet(correlated_fitness_surface(prep, "survival", traits, method = "gam", grid_n = 20,
                                               mask = TRUE, group = "sex", group_effect = TRUE, k = NULL,
                                               bs = "tp", smoothing = "REML", clamp = TRUE,
                                               count_family = "poisson"))
      fit <<- list(prep = prep, model = surf$model, traits = traits)
    }
    fit
  }
})

# the simulate loop exactly as it was before the methods were added
reference_landscape <- function(data, model, traits, Sigma, grid, extra) {
  vapply(seq_len(nrow(grid)), function(i) {
    simulated <- MASS::mvrnorm(n = 50, mu = as.numeric(grid[i, ]), Sigma = Sigma)
    colnames(simulated) <- traits
    sim_df <- as.data.frame(simulated)
    for (v in names(extra)) sim_df[[v]] <- extra[[v]]
    mean(predict(model, newdata = sim_df, type = "response"), na.rm = TRUE)
  }, numeric(1))
}

test_that("simulate gives the numbers it gave before the methods were added", {
  b <- bumpus_surface()
  set.seed(101)
  land <- quiet(adaptive_landscape(b$prep, b$model, b$traits, group_col = "sex", grid_n = 10, simulation_n = 50))
  after <- runif(1)
  expect_identical(land$method, "simulate")
  expect_null(land$lattice)
  expect_equal(land$simulation_n, 50)

  # the same loop written out, from the same seed
  set.seed(101)
  ref <- reference_landscape(b$prep, b$model, b$traits, land$population_variance,
                             land$grid[, b$traits], list(sex = "male"))
  expect_identical(land$grid$.mean_fit, ref)
  expect_identical(runif(1), after)

  # and the values the old loop gave on the machine the methods were written on;
  # the draws depend on the platform's eigen decomposition, so only there
  skip_on_cran()
  skip_if_not(Sys.info()[["sysname"]] == "Darwin" && R.version$arch == "aarch64")
  expect_equal(land$grid$.mean_fit[c(1, 17, 55, 100)],
               c(0.9423220666309048, 0.7240915544760781, 0.5336278195436012, 0.0690873253861986),
               tolerance = 1e-8)
  expect_equal(sum(land$grid$.mean_fit), 53.3231147316238, tolerance = 1e-8)
  expect_equal(land$grid$.outside[c(1, 17, 55, 100)], c(0.92, 1, 0.06, 0.96))
  expect_equal(after, 0.0561556394677609, tolerance = 1e-12)
})

test_that("the lattice follows simulate draw for draw", {
  b <- bumpus_surface()
  set.seed(7)
  sim <- quiet(adaptive_landscape(b$prep, b$model, b$traits, group_col = "sex", grid_n = 12, simulation_n = 200))
  rng_sim <- .Random.seed
  set.seed(7)
  msgs <- capture_messages(
    lat <- suppressWarnings(adaptive_landscape(b$prep, b$model, b$traits, group_col = "sex", grid_n = 12,
                                               simulation_n = 200, method = "lattice", lattice_n = 100))
  )
  expect_match(msgs, "largest error .* within the tolerance", all = FALSE)
  expect_identical(.Random.seed, rng_sim)
  expect_identical(lat$method, "lattice")
  expect_lt(max(abs(lat$grid$.mean_fit - sim$grid$.mean_fit)), 1e-4)
  expect_identical(rownames(lat$optimum), rownames(sim$optimum))
  expect_equal(lat$grid$.outside, sim$grid$.outside)
  expect_equal(lat$support, sim$support)
  expect_true(lat$lattice_check$passed)
  expect_lte(lat$lattice_check$n_checked, 2000)
  expect_gt(lat$lattice_check$n_checked, 1500)
  expect_equal(unname(lat$lattice$n), c(100, 100))
  expect_output(print(lat), "Method: lattice")
  expect_output(print(lat), "Lattice check")
})

test_that("a coarse lattice is made finer, and the draws predicted directly when that costs less", {
  set.seed(2)
  df <- data.frame(z1 = as.numeric(scale(rnorm(300))), z2 = as.numeric(scale(rnorm(300))))
  df$w <- sin(df$z1) + cos(df$z2) + rnorm(300, 0, 0.05)
  df$v <- sin(2 * df$z1) + cos(2 * df$z2) + rnorm(300, 0, 0.05)
  fit <- mgcv::gam(w ~ s(z1, z2, k = 40), data = df)
  sharp <- mgcv::gam(v ~ s(z1, z2, k = 40), data = df)

  # 10 x 10 nodes are far too coarse for this surface; 160 x 160 pass
  set.seed(3)
  sim <- quiet(adaptive_landscape(df, fit, c("z1", "z2"), grid_n = 10, simulation_n = 500))
  set.seed(3)
  fine <- quiet(adaptive_landscape(df, fit, c("z1", "z2"), grid_n = 10, simulation_n = 500,
                                   method = "lattice", lattice_n = 10))
  expect_gt(fine$lattice_check$refinements, 0)
  expect_true(fine$lattice_check$passed)
  expect_equal(fine$lattice_check$direct, 0L)
  expect_gt(fine$lattice_check$lattice_n, 10)
  expect_lt(max(abs(fine$grid$.mean_fit - sim$grid$.mean_fit)), 1e-3 * diff(range(sim$grid$.mean_fit)))

  # for the sharper surface 2500 nodes fail the check, and 10000 would be
  # more than the 5000 draws
  set.seed(3)
  small <- quiet(adaptive_landscape(df, sharp, c("z1", "z2"), grid_n = 10, simulation_n = 50))
  set.seed(3)
  msgs <- capture_messages(
    direct <- adaptive_landscape(df, sharp, c("z1", "z2"), grid_n = 10, simulation_n = 50,
                                 method = "lattice", lattice_n = 50)
  )
  expect_match(msgs, "predicted directly", all = FALSE)
  expect_equal(direct$lattice_check$direct, 5000L)
  expect_equal(direct$grid$.mean_fit, small$grid$.mean_fit, tolerance = 1e-10)
  expect_output(print(direct), "5000 draws predicted directly")

  # 400 nodes are far off, and the 1600 + 6400 a pass should take would cost
  # more than the 5000 draws: they are predicted directly at once, rather
  # than after a refinement that cannot pay
  set.seed(3)
  at_once <- quiet(adaptive_landscape(df, sharp, c("z1", "z2"), grid_n = 10, simulation_n = 50,
                                      method = "lattice", lattice_n = 20))
  expect_equal(at_once$lattice_check$refinements, 0L)
  expect_equal(at_once$lattice$predictions, 400)
  expect_equal(at_once$lattice_check$direct, 5000L)
  expect_equal(at_once$grid$.mean_fit, small$grid$.mean_fit, tolerance = 1e-10)
})

test_that("the check's tolerance follows fitness over the grid, not the lattice's margin", {
  # count fitness under directional selection: the log link climbs far beyond
  # the data, so fitness out in the lattice's margin dwarfs the landscape
  set.seed(11)
  df <- data.frame(z1 = as.numeric(scale(rnorm(300))), z2 = as.numeric(scale(rnorm(300))))
  df$kids <- rpois(300, exp(0.2 + 0.9 * df$z1 + 0.6 * sin(1.5 * df$z2)))
  fit <- mgcv::gam(kids ~ s(z1, z2), family = poisson, data = df, method = "REML")
  set.seed(1)
  sim <- quiet(adaptive_landscape(df, fit, c("z1", "z2"), grid_n = 12, simulation_n = 200))
  set.seed(1)
  lat <- quiet(adaptive_landscape(df, fit, c("z1", "z2"), grid_n = 12, simulation_n = 200, method = "lattice",
                                  lattice_n = 8))
  span <- diff(range(sim$grid$.mean_fit))
  # a thousandth of fitness over the grid is about a thousandth of the landscape
  expect_lt(lat$lattice_check$tol, 2e-3 * span)
  # so 8 x 8 nodes are refined, and the landscape lands within 0.1% of simulate
  expect_gt(lat$lattice_check$refinements, 0)
  expect_true(lat$lattice_check$passed)
  expect_lt(max(abs(lat$grid$.mean_fit - sim$grid$.mean_fit)), 1e-3 * span)
})

test_that("lattice_n is checked only where it is used", {
  b <- bumpus_surface()
  expect_error(adaptive_landscape(b$prep, b$model, b$traits, group_col = "sex", grid_n = 5, simulation_n = 10,
                                  method = "lattice", lattice_n = 4), "lattice_n")
  expect_no_warning(suppressMessages(adaptive_landscape(b$prep, b$model, b$traits, group_col = "sex", grid_n = 5,
                                                        simulation_n = 10, lattice_n = Inf)))
})

test_that("the lattice falls back to simulate when it would cost more than the draws", {
  b <- bumpus_surface()
  set.seed(9)
  sim <- quiet(adaptive_landscape(b$prep, b$model, b$traits, group_col = "sex", grid_n = 5, simulation_n = 20))
  set.seed(9)
  msgs <- capture_messages(
    lat <- adaptive_landscape(b$prep, b$model, b$traits, group_col = "sex", grid_n = 5, simulation_n = 20,
                              method = "lattice")
  )
  expect_match(msgs, "predicted directly", all = FALSE)
  expect_identical(lat$method, "simulate")
  expect_identical(lat$grid, sim$grid)
})

test_that("convolve matches the closed form for a quadratic fitness function", {
  # w(z) = a + b'z + z'Gz/2 averaged over N(mu, Sigma) is
  # a + b'mu + (mu'G mu + tr(G Sigma))/2
  set.seed(11)
  z <- MASS::mvrnorm(300, c(0, 0), matrix(c(1, 0.6, 0.6, 1), 2))
  df <- data.frame(z1 = as.numeric(scale(z[, 1])), z2 = as.numeric(scale(z[, 2])))
  a <- 1
  b <- c(0.3, -0.2)
  G <- matrix(c(-0.4, 0.15, 0.15, -0.3), 2)
  Z <- as.matrix(df)
  df$y <- as.numeric(a + Z %*% b + 0.5 * rowSums((Z %*% G) * Z))
  # a thin-plate smooth with m = 3 leaves every quadratic unpenalised, so it
  # reproduces this one exactly, also beyond the data
  fit <- mgcv::gam(y ~ s(z1, z2, m = 3), data = df)

  land <- quiet(adaptive_landscape(df, fit, c("z1", "z2"), grid_n = 9, method = "convolve"))
  S <- land$population_variance
  expect_gt(abs(S[1, 2]), 0.4)
  mu <- as.matrix(land$grid[, c("z1", "z2")])
  exact <- as.numeric(a + mu %*% b + 0.5 * (rowSums((mu %*% G) * mu) + sum(diag(G %*% S))))
  expect_lt(max(abs(land$grid$.mean_fit - exact)), 1e-5)
  expect_identical(land$method, "convolve")
  expect_true(is.na(land$simulation_n))
  expect_null(land$lattice_check)
  # every grid point is a lattice node, at most a quarter SD apart
  expect_true(all(land$lattice$spacing <= sqrt(diag(S)) / 4 + 1e-12))
  expect_output(print(land), "Method: convolve")

  # cubic interpolation reproduces a quadratic, so the lattice method gives
  # what predicting every draw gives
  set.seed(5)
  sim <- quiet(adaptive_landscape(df, fit, c("z1", "z2"), grid_n = 9, simulation_n = 400))
  set.seed(5)
  lat <- quiet(adaptive_landscape(df, fit, c("z1", "z2"), grid_n = 9, simulation_n = 400, method = "lattice",
                                  lattice_n = 50))
  expect_lt(lat$lattice_check$max_error, 1e-8)
  expect_lt(max(abs(lat$grid$.mean_fit - sim$grid$.mean_fit)), 1e-8)

  # nearly collinear traits squeeze the distribution into a narrow ridge; a
  # lattice a quarter SD apart along each trait missed it by 1e-4 at r = 0.99
  # and 6e-4 at r = 0.995, so the spacing follows the conditional SD
  for (r in c(0.99, -0.995)) {
    S_r <- matrix(c(1, r, r, 1), 2)
    near <- quiet(adaptive_landscape(df, fit, c("z1", "z2"), population_variance = S_r, grid_n = 17,
                                     custom_range = list(z1 = c(-2, 2), z2 = c(-2, 2)), method = "convolve"))
    expect_identical(near$method, "convolve")
    mu <- as.matrix(near$grid[, c("z1", "z2")])
    exact <- as.numeric(a + mu %*% b + 0.5 * (rowSums((mu %*% G) * mu) + sum(diag(G %*% S_r))))
    expect_lt(max(abs(near$grid$.mean_fit - exact)), 2e-5)
  }
})

test_that("convolve hands over to the lattice when its lattice would cost too much", {
  b <- bumpus_surface()
  # a variance that is not positive definite becomes a singular one, all its
  # weight on a line no rectangular lattice can follow; convolve used to
  # return the individual fitness surface itself there
  bad <- matrix(c(1, 1.2, 1.2, 1), 2)
  set.seed(6)
  msgs <- capture_messages(
    conv <- adaptive_landscape(b$prep, b$model, b$traits, group_col = "sex", grid_n = 6, simulation_n = 100,
                               population_variance = bad, method = "convolve", lattice_n = 20)
  )
  expect_match(msgs, "using the lattice method instead", all = FALSE)
  expect_identical(conv$method, "lattice")
  set.seed(6)
  lat <- quiet(adaptive_landscape(b$prep, b$model, b$traits, group_col = "sex", grid_n = 6, simulation_n = 100,
                                  population_variance = bad, method = "lattice", lattice_n = 20))
  expect_identical(conv$grid, lat$grid)
  w_mu <- as.numeric(predict(b$model, data.frame(conv$grid[, b$traits], sex = "male"), type = "response"))
  expect_gt(max(abs(conv$grid$.mean_fit - w_mu)), 0.05)

  # a tiny within-population SD next to the grid step: hundreds of thousands
  # of nodes for 2000 draws
  set.seed(6)
  msgs <- capture_messages(
    tiny <- adaptive_landscape(b$prep, b$model, b$traits, group_col = "sex", grid_n = 20, simulation_n = 100,
                               population_variance = diag(c(1e-4, 1e-4)), method = "convolve", lattice_n = 20)
  )
  expect_match(msgs, "using the lattice method instead", all = FALSE)
  expect_identical(tiny$method, "lattice")
})

test_that("convolve reaches far enough for fitness that climbs steeply under a log link", {
  # E[exp(a + b z)] for z ~ N(mu, 1) is exp(a + b mu + b^2 / 2); with b = 2.5
  # the integrand peaks 2.5 SD out, and weights cut at 5 SD lost 0.25% of it
  d <- data.frame(z = as.numeric(scale(seq(-2, 2, length.out = 200))))
  d$y <- exp(0.2 + 2.5 * d$z)
  fit <- suppressWarnings(mgcv::gam(y ~ s(z), family = quasipoisson, data = d))
  conv <- quiet(adaptive_landscape(d, fit, "z", population_variance = matrix(1, 1, 1, dimnames = list("z", "z")),
                                   grid_n = 11, method = "convolve"))
  exact <- exp(0.2 + 2.5 * conv$grid$z + 2.5^2 / 2)
  expect_lt(max(abs(conv$grid$.mean_fit / exact - 1)), 1e-6)
  expect_equal(conv$lattice$reach_sd, 8)
})

test_that("convolve agrees with a large simulation", {
  b <- bumpus_surface()
  conv <- quiet(adaptive_landscape(b$prep, b$model, b$traits, group_col = "sex", grid_n = 5, method = "convolve"))
  # the lattice method follows simulate draw for draw (tested above) and is far quicker
  set.seed(4)
  sim <- quiet(adaptive_landscape(b$prep, b$model, b$traits, group_col = "sex", grid_n = 5, simulation_n = 20000,
                                  method = "lattice"))
  # survival lies in [0, 1], so the Monte Carlo SE is at most sqrt(W(1 - W) / n)
  se <- sqrt(conv$grid$.mean_fit * (1 - conv$grid$.mean_fit) / 20000)
  expect_true(all(abs(sim$grid$.mean_fit - conv$grid$.mean_fit) < 4 * se + 1e-4))
  # a share has SE at most 0.5 / sqrt(20000), 0.0035
  expect_lt(max(abs(sim$grid$.outside - conv$grid$.outside)), 0.015)
})

test_that("convolve does not depend on the seed and leaves the random numbers alone", {
  b <- bumpus_surface()
  set.seed(1)
  msgs <- capture_messages(
    one <- adaptive_landscape(b$prep, b$model, b$traits, group_col = "sex", grid_n = 8, method = "convolve")
  )
  rng <- .Random.seed
  # simulation_n is only mentioned when it was given
  expect_false(any(grepl("simulation_n", msgs)))
  set.seed(2)
  msgs <- capture_messages(
    two <- adaptive_landscape(b$prep, b$model, b$traits, group_col = "sex", grid_n = 8, simulation_n = 10,
                              method = "convolve")
  )
  expect_match(msgs, "simulation_n is not used", all = FALSE)
  expect_identical(one$grid, two$grid)
  expect_identical(one$optimum, two$optimum)
  set.seed(1)
  quiet(adaptive_landscape(b$prep, b$model, b$traits, group_col = "sex", grid_n = 8, method = "convolve"))
  expect_identical(.Random.seed, rng)
})

test_that("both methods draw one-trait landscapes", {
  b <- bumpus_surface()
  uni <- quiet(univariate_spline(b$prep, "survival", "weight", fitness_type = "binary", group = "sex", k = 5,
                                 bs = "cr", smoothing = "GCV.Cp", bootstrap = FALSE))
  set.seed(5)
  sim <- quiet(adaptive_landscape(b$prep, uni$model, "weight", group_col = "sex", grid_n = 20, simulation_n = 200))
  rng_sim <- .Random.seed
  set.seed(5)
  lat <- quiet(adaptive_landscape(b$prep, uni$model, "weight", group_col = "sex", grid_n = 20, simulation_n = 200,
                                  method = "lattice", lattice_n = 50))
  expect_identical(.Random.seed, rng_sim)
  expect_equal(unname(lat$lattice$n), 1000)
  expect_lt(max(abs(lat$grid$.mean_fit - sim$grid$.mean_fit)), 1e-4)
  expect_equal(lat$grid$.outside, sim$grid$.outside)
  expect_equal(lat$grid$.ind_fit, sim$grid$.ind_fit)

  conv <- quiet(adaptive_landscape(b$prep, uni$model, "weight", group_col = "sex", grid_n = 20, method = "convolve"))
  expect_equal(conv$grid$.ind_fit, sim$grid$.ind_fit)
  # against the integral done numerically
  sdv <- sqrt(as.matrix(conv$population_variance)[1, 1])
  w <- function(x) as.numeric(predict(uni$model, data.frame(weight = x, sex = "male"), type = "response"))
  exact <- vapply(conv$grid$weight, function(m) {
    stats::integrate(function(x) w(x) * stats::dnorm(x, m, sdv), m - 8 * sdv, m + 8 * sdv, rel.tol = 1e-10)$value
  }, numeric(1))
  expect_lt(max(abs(conv$grid$.mean_fit - exact)), 1e-6)
  # the share of the population beyond the observed range, from the part of
  # each lattice cell beyond it (taking whole nodes was off by up to 0.006)
  r <- range(b$prep$weight)
  beyond <- stats::pnorm(r[1], conv$grid$weight, sdv) + stats::pnorm(r[2], conv$grid$weight, sdv, lower.tail = FALSE)
  expect_lt(max(abs(conv$grid$.outside - beyond)), 1e-3)
})

test_that("both methods clamp a thin-plate surface", {
  skip_if_not_installed("fields")
  set.seed(3)
  df <- data.frame(z1 = as.numeric(scale(rnorm(120))), z2 = as.numeric(scale(rnorm(120))))
  df$w <- rbinom(120, 1, plogis(2 * df$z1))
  held <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 12, method = "tps", mask = FALSE))

  set.seed(8)
  sim <- quiet(adaptive_landscape(df, held$model, c("z1", "z2"), simulation_n = 300, grid_n = 8))
  set.seed(8)
  lat <- quiet(adaptive_landscape(df, held$model, c("z1", "z2"), simulation_n = 300, grid_n = 8,
                                  method = "lattice", lattice_n = 60))
  expect_true(lat$clamp)
  expect_gt(lat$clipped, 0)
  expect_equal(lat$clipped, sim$clipped)
  expect_lt(max(abs(lat$grid$.mean_fit - sim$grid$.mean_fit)), 1e-4)
  expect_true(all(lat$grid$.mean_fit >= 0 & lat$grid$.mean_fit <= 1))

  conv <- quiet(adaptive_landscape(df, held$model, c("z1", "z2"), grid_n = 4, method = "convolve"))
  expect_gt(conv$clipped, 0)
  expect_true(all(conv$grid$.mean_fit >= 0 & conv$grid$.mean_fit <= 1))
  set.seed(9)
  big <- quiet(adaptive_landscape(df, held$model, c("z1", "z2"), simulation_n = 20000, grid_n = 4,
                                  method = "lattice"))
  se <- sqrt(conv$grid$.mean_fit * (1 - conv$grid$.mean_fit) / 20000)
  expect_true(all(abs(big$grid$.mean_fit - conv$grid$.mean_fit) < 4 * se + 1e-3))
  # without clamping the surface leaves 0 to 1 and so does its average
  raw <- quiet(adaptive_landscape(df, held$model, c("z1", "z2"), grid_n = 4, method = "convolve", clamp = FALSE))
  expect_equal(raw$clipped, 0L)
  expect_gte(max(raw$grid$.mean_fit), max(conv$grid$.mean_fit))
})
