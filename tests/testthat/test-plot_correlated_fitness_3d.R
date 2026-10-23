# draw to a null device so the tests leave no files behind
draws <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expr
}

test_that("the 3D surface needs two traits", {
  cfs <- list(
    grid = data.frame(z1 = rep(1:5, 5), z2 = rep(1:5, each = 5), .fit = runif(25))
  )
  expect_error(plot_correlated_fitness_3d(cfs, "z1"), "needs two traits")
  expect_error(plot_correlated_fitness_3d(cfs, c("z1", "z3")), "not found")
})

test_that("an incomplete grid is refused", {
  cfs <- list(
    grid = data.frame(z1 = rep(1:5, 5), z2 = rep(1:5, each = 5), .fit = runif(25))[-3, ]
  )
  expect_error(plot_correlated_fitness_3d(cfs, c("z1", "z2")), "not a complete regular grid")
})

test_that("tps and gam surfaces draw in 3D, masked cells left out", {
  set.seed(3)
  n <- 200
  df <- data.frame(z1 = rnorm(n), z2 = rnorm(n))
  df$w <- 1 - 0.3 * df$z1^2 - 0.2 * df$z2^2 + rnorm(n, 0, 0.1)
  tps <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 20, method = "tps"))
  gam <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 20, method = "gam"))
  expect_true(anyNA(tps$grid$.fit))
  pm <- draws(plot_correlated_fitness_3d(tps, c("z1", "z2")))
  expect_equal(dim(pm), c(4L, 4L))
  expect_no_error(draws(plot_correlated_fitness_3d(gam, c("z1", "z2"), theta = 20, phi = 15,
                                                   color_palette = grDevices::grey.colors(20))))
  expect_no_error(draws(plot_correlated_fitness_3d(gam, c("z1", "z2"), show_peaks = FALSE)))
})

test_that("the 3D surface works with groups", {
  set.seed(6)
  n <- 240
  df <- data.frame(g = rep(c("a", "b"), each = n / 2), z1 = rnorm(n), z2 = rnorm(n))
  df$w <- rbinom(n, 1, plogis(1 - 0.8 * df$z1^2))
  prep <- quiet(prepare_selection_data(df, "w", c("z1", "z2"), group = "g"))
  one <- quiet(correlated_fitness_surface(prep, "w", c("z1", "z2"), grid_n = 15,
                                          group = "g", group_effect = FALSE))
  expect_false(is.null(one$groups))
  expect_no_error(draws(plot_correlated_fitness_3d(one, c("z1", "z2"))))
  each <- quiet(correlated_fitness_surface(prep, "w", c("z1", "z2"), grid_n = 15,
                                           group = "g", by_group = TRUE))
  for (s in each) expect_no_error(draws(plot_correlated_fitness_3d(s, c("z1", "z2"))))
})
