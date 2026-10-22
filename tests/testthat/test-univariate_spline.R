library(testthat)

test_that("univariate_spline generates valid splines", {
  set.seed(42)
  df <- data.frame(
    w = rnorm(50, 1, 0.1),
    z = rnorm(50)
  )
  
  df$z <- as.numeric(scale(df$z))
  
  res <- suppressWarnings(univariate_spline(df, "w", "z", fitness_type = "continuous", k = 3))
  expect_s3_class(res, "univariate_fitness")
  expect_true(all(c("fit", "lwr", "upr") %in% names(res$grid)))
})

test_that("the basis and smoothing method can be chosen and the defaults are unchanged", {
  set.seed(6)
  df <- data.frame(z = as.numeric(scale(rnorm(150))))
  df$w <- 1 + 0.3 * df$z - 0.25 * df$z^2 + rnorm(150, 0, 0.2)

  default <- suppressMessages(univariate_spline(df, "w", "z", fitness_type = "continuous", k = 6))
  explicit <- suppressMessages(univariate_spline(df, "w", "z", fitness_type = "continuous", k = 6, bs = "cr", smoothing = "GCV.Cp"))
  expect_equal(default$grid$fit, explicit$grid$fit)
  expect_equal(default$basis, "cr")
  expect_equal(default$smoothing, "GCV.Cp")
  expect_equal(default$spline_type, "cubic regression spline (GCV)")

  tp <- suppressMessages(univariate_spline(df, "w", "z", fitness_type = "continuous", k = 6, bs = "tp", smoothing = "REML"))
  expect_equal(tp$basis, "tp")
  expect_equal(tp$model$method, "REML")
  expect_equal(tp$spline_type, "thin-plate spline (REML)")
  # a different smoother, but the same shape of fitness function
  expect_gt(cor(default$grid$fit, tp$grid$fit), 0.99)

  expect_error(univariate_spline(df, "w", "z", fitness_type = "continuous", bs = "ad"), "should be one of")
})

test_that("count fitness is fitted on the raw counts with a Poisson family", {
  set.seed(5)
  df <- data.frame(z = as.numeric(scale(rnorm(120))))
  df$kids <- rpois(120, exp(0.3 + 0.4 * df$z))

  res <- suppressMessages(univariate_spline(df, "kids", "z", k = 4))
  expect_equal(res$fitness_type, "count")
  expect_equal(res$family, "poisson(log)")
  expect_true(all(res$grid$fit >= 0))
  expect_true(all(res$grid$lwr >= 0))
})

test_that("refits on a cluster give the same ribbon and random number state as the serial loop", {
  skip_on_cran()
  skip_if_not_installed("parallel")
  set.seed(8)
  df <- data.frame(z = as.numeric(scale(rnorm(80))))
  df$s <- rbinom(80, 1, plogis(0.4 - 0.5 * df$z^2))
  df$kids <- rpois(80, exp(0.3 + 0.2 * df$z))
  cl <- parallel::makePSOCKcluster(2)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  run <- function(fitness, ...) {
    set.seed(3)
    fit <- quiet(univariate_spline(df, fitness, "z", k = 5, bootstrap = TRUE, n_boot = 20, ...))
    list(grid = fit$grid, ci = fit$ci_method, seed = .Random.seed)
  }
  serial <- run("s")
  expect_identical(serial$ci, "bootstrap (percentile)")
  expect_identical(run("s", cl = cl), serial)
  # a negative binomial fit refits serially, cluster or not
  expect_identical(run("kids", count_family = "nb", cl = cl), run("kids", count_family = "nb"))
  expect_error(univariate_spline(df, "s", "z", bootstrap = TRUE, cl = 2), "cluster")
  # the workers fit under this session's contrasts
  df$g <- rep(c("a", "b", "c"), length.out = 80)
  op <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(op), add = TRUE)
  expect_identical(run("kids", group = "g", cl = cl), run("kids", group = "g"))
  options(op)
  # with options(warn = 2) a refit's warning is an error that drops the refit,
  # with or without workers (few survivors, so some refits warn)
  set.seed(1)
  rare <- data.frame(z = as.numeric(scale(rnorm(40))))
  rare$s <- as.integer(seq_len(40) %in% sample.int(40, 4))
  strict <- function(...) {
    old <- options(warn = 2)
    on.exit(options(old))
    set.seed(1)
    out <- tryCatch(suppressMessages(univariate_spline(rare, "s", "z", k = 5, bootstrap = TRUE, n_boot = 40, ...)),
                    error = conditionMessage)
    list(out = if (is.character(out)) out else out$grid, seed = .Random.seed)
  }
  expect_identical(strict(cl = cl), strict())
})