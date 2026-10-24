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

test_that("each grid cell is outlined, and NA leaves the outlines off", {
  skip_if_not_installed("fields")
  set.seed(4)
  n <- 200
  df <- data.frame(z1 = rnorm(n), z2 = rnorm(n))
  df$w <- 1 - 0.3 * df$z1^2 - 0.2 * df$z2^2 + rnorm(n, 0, 0.1)
  s <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 20, method = "tps"))
  seen <- NULL
  local_mocked_bindings(drape.plot = function(x, y, z, ...) {
    seen <<- list(n = c(length(x), length(y)), z = z, ...)
    diag(4)
  }, .package = "fields")
  plot_correlated_fitness_3d(s, c("z1", "z2"), show_peaks = FALSE)
  expect_equal(seen$border, "grey30")
  expect_equal(seen$lwd, 0.5)
  expect_equal(seen$n, c(20L, 20L))
  # masked cells stay NA, so persp() draws neither the cell nor its outline
  expect_equal(sum(is.na(seen$z)), sum(is.na(s$grid$.fit)))
  plot_correlated_fitness_3d(s, c("z1", "z2"), show_peaks = FALSE, border = NA)
  expect_true(is.na(seen$border))
  plot_correlated_fitness_3d(s, c("z1", "z2"), show_peaks = FALSE, border = "grey40", lwd = 1)
  expect_equal(seen$border, "grey40")
  expect_equal(seen$lwd, 1)
})

test_that("the shading is light by default and shade passes through", {
  skip_if_not_installed("fields")
  set.seed(6)
  n <- 200
  df <- data.frame(z1 = rnorm(n), z2 = rnorm(n))
  df$w <- 1 - 0.3 * df$z1^2 - 0.2 * df$z2^2 + rnorm(n, 0, 0.1)
  s <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 20, method = "tps"))
  seen <- NULL
  local_mocked_bindings(drape.plot = function(x, y, z, ...) {
    seen <<- list(...)
    diag(4)
  }, .package = "fields")
  plot_correlated_fitness_3d(s, c("z1", "z2"), show_peaks = FALSE)
  expect_equal(seen$shade, 0.15)
  plot_correlated_fitness_3d(s, c("z1", "z2"), show_peaks = FALSE, shade = 0.5)
  expect_equal(seen$shade, 0.5)
  plot_correlated_fitness_3d(s, c("z1", "z2"), show_peaks = FALSE, shade = NA)
  expect_true(is.na(seen$shade))
})

test_that("the surface draws with shade off, light and dark", {
  skip_if_not_installed("fields")
  set.seed(7)
  n <- 200
  df <- data.frame(z1 = rnorm(n), z2 = rnorm(n))
  df$w <- 1 - 0.3 * df$z1^2 - 0.2 * df$z2^2 + rnorm(n, 0, 0.1)
  s <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 20, method = "gam"))
  for (sh in list(NA, 0, 0.15, 0.5)) {
    expect_no_error(draws(plot_correlated_fitness_3d(s, c("z1", "z2"), shade = sh)))
  }
})

test_that("the default outlines are dark grey, fading only on very fine grids", {
  expect_equal(formals(plot_correlated_fitness_3d)$border, "grey30")
  expect_equal(grid_line_colour("grey30", 80), "grey30")
  expect_equal(grid_line_colour("black", 20), "black")
  expect_equal(grid_line_colour("black", 100), "black")
  expect_true(is.na(grid_line_colour(NA, 200)))
  fine <- grDevices::col2rgb(grid_line_colour("black", 200), alpha = TRUE)
  finer <- grDevices::col2rgb(grid_line_colour("black", 400), alpha = TRUE)
  expect_lt(fine["alpha", 1], 255)
  expect_lt(finer["alpha", 1], fine["alpha", 1])
  grey <- grDevices::col2rgb(grid_line_colour("grey30", 200), alpha = TRUE)
  expect_equal(unname(grey[1:3, 1]), c(77, 77, 77))
  expect_lt(grey["alpha", 1], 255)
})

test_that("an 80 by 80 surface draws with its grid lines, fine grids too", {
  skip_if_not_installed("fields")
  set.seed(5)
  n <- 200
  df <- data.frame(z1 = rnorm(n), z2 = rnorm(n))
  df$w <- 1 - 0.3 * df$z1^2 - 0.2 * df$z2^2 + rnorm(n, 0, 0.1)
  for (k in c(5, 80, 120)) {
    s <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = k, method = "gam"))
    expect_no_error(draws(plot_correlated_fitness_3d(s, c("z1", "z2"))))
  }
  s <- quiet(correlated_fitness_surface(df, "w", c("z1", "z2"), grid_n = 120, method = "gam"))
  expect_no_error(draws(plot_correlated_fitness_3d(s, c("z1", "z2"), border = "black")))
})
