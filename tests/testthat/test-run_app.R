library(testthat)

test_that("the app ships with the package and run_app points at it", {
  dir <- system.file("app", package = "lande")
  expect_true(nzchar(dir))
  expect_true(file.exists(file.path(dir, "app.R")))
  expect_true(is.function(run_app))
  # the app file parses and refers to the package
  expr <- parse(file.path(dir, "app.R"))
  expect_gt(length(expr), 10)
  calls <- vapply(expr, function(e) paste(deparse(e), collapse = ""), "")
  expect_true(any(grepl("library(lande)", calls, fixed = TRUE)))
})

test_that("the app's jet theme colours the 2D bands and the 3D surface", {
  dir <- system.file("app", package = "lande")
  expr <- parse(file.path(dir, "app.R"))
  # only the colour definitions, without starting the app
  keep <- vapply(expr, function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                   deparse(e[[2]]) %in% c("jet", "THEMES"), logical(1))
  env <- new.env()
  for (e in expr[keep]) eval(e, env)
  th <- env$THEMES[["Rainbow (jet), white and black"]]
  expect_false(is.null(th))
  expect_equal(th$ramp(9)[c(1, 5, 9)], c("#00007F", "#7FFF7F", "#7F0000"))
  expect_equal(th$ramp(100), jet_colours(100))
  expect_equal(th$points, c("white", "black"))
  d <- expand.grid(x = 1:20, y = 1:20)
  d$z <- with(d, -(x - 10)^2 - (y - 10)^2)
  p <- ggplot2::ggplot(d, ggplot2::aes(x, y, z = z)) + ggplot2::geom_contour_filled(bins = 12) +
    th$fill(name = "Fitness")
  # every band is coloured, from dark blue at the bottom to dark red at the top
  ld <- ggplot2::layer_data(p)
  n <- length(unique(ld$level))
  expect_false(anyNA(ld$fill))
  expect_equal(ld$fill[order(ld$level)][c(1, nrow(ld))], env$jet(n)[c(1, n)])
})
