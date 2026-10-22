# Regression tests for behaviour that broke, or nearly broke, during the
# package restructure. Each block names the failure it guards against.

sim_data <- function(n = 200, seed = 11) {
  set.seed(seed)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  data.frame(
    w = 5 + 1.5 * z1 + 0.4 * z1^2 + rnorm(n, 0, 0.5),
    z1 = z1, z2 = z2,
    grp = rep(c("A", "B"), length.out = n)
  )
}

test_that("direct analyzer calls fit the column they are given", {
  df <- sim_data()
  prep <- suppressWarnings(suppressMessages(prepare_selection_data(df, "w", c("z1", "z2"))))

  ref <- suppressWarnings(suppressMessages(
    selection_coefficients(df, "w", c("z1", "z2"), fitness_type = "continuous")))
  lin <- suppressWarnings(suppressMessages(
    analyze_linear_selection(prep, "relative_fitness", c("z1", "z2"), "continuous")))
  nl <- suppressWarnings(suppressMessages(
    analyze_nonlinear_selection(prep, "relative_fitness", c("z1", "z2"), "continuous")))
  # the raw column is fitted as asked, with a note about the relative one
  expect_message(
    raw <- suppressWarnings(analyze_linear_selection(prep, "w", c("z1", "z2"), "continuous")),
    "relative_fitness"
  )
  expect_equal(unname(coef(raw$model)["z1"]), unname(coef(lm(w ~ z1 + z2, data = prep))["z1"]))

  expect_equal(unname(coef(lin$model)["z1"]),
               ref$Beta_Coefficient[ref$Term == "z1"], tolerance = 1e-10)
  expect_equal(2 * unname(coef(nl$model)["I(z1^2)"]),
               ref$Beta_Coefficient[ref$Term == "z1²"], tolerance = 1e-10)
})

test_that("a stale relative_fitness column does not change the gradients", {
  set.seed(7)
  d <- data.frame(g = rep(c("a", "b"), each = 100), z1 = rnorm(200), z2 = rnorm(200))
  d$w <- ifelse(d$g == "a", 2, 6) + 0.5 * d$z1 + rnorm(200, 0, 0.5)
  # the same data with a leftover relative_fitness column, worked out within groups
  stale <- d
  stale$relative_fitness <- d$w / ave(d$w, d$g)
  fit <- function(x, ...) suppressWarnings(suppressMessages(selection_coefficients(x, "w", c("z1", "z2"), ...)))
  expect_equal(fit(stale)$Beta_Coefficient, fit(d)$Beta_Coefficient, tolerance = 1e-8)
  expect_equal(fit(stale, use_relative_for_fit = FALSE)$Beta_Coefficient,
               fit(d, use_relative_for_fit = FALSE)$Beta_Coefficient, tolerance = 1e-8)
})

test_that("a direct binary call with a non-0/1 column is rejected clearly", {
  d <- data.frame(surv = rbinom(60, 1, 0.5), z = rnorm(60))
  prep <- suppressWarnings(suppressMessages(prepare_selection_data(d, "surv", "z")))
  expect_error(
    analyze_linear_selection(prep, "relative_fitness", "z", "binary"),
    "raw 0/1 column"
  )
})

test_that("single-trait analyses do not emit a VIF warning", {
  df <- sim_data()
  msgs <- character(0)
  withCallingHandlers(
    analyze_disruptive_selection(df, "w", "z1", fitness_type = "continuous"),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage")
  )
  expect_false(any(grepl("VIF", msgs)))
})

test_that("an exact least-squares fit does not stop selection_coefficients()", {
  # 3 birds for the 3 coefficients of the quadratic model: no residual degrees
  # of freedom, so its VIFs are NaN, which the collinearity check used to trip on
  skip_if_not_installed("car")
  d <- data.frame(w = c(1, 2, 0.5), x = c(0.1, 0.7, 1.9))
  msgs <- character(0)
  res <- withCallingHandlers(
    suppressMessages(selection_coefficients(d, "w", "x", fitness_type = "continuous")),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(res$Term, c("x", "x\u00b2"))
  expect_true(all(is.finite(res$Beta_Coefficient)))
  expect_true(is.nan(res$Standard_Error[2]))
  expect_false(any(grepl("Collinear", msgs)))
})

test_that("rows with a missing group label are kept and standardised", {
  df <- sim_data(n = 40)
  df$grp[c(3, 7, 11)] <- NA
  out <- suppressWarnings(suppressMessages(
    prepare_selection_data(df, "w", c("z1", "z2"), group = "grp", na_action = "warn")))
  expect_equal(nrow(out), 40)
  expect_false(anyNA(out$z1))
  expect_false(anyNA(out$relative_fitness))
  # the labelled groups are still standardised on their own
  expect_equal(sd(out$z1[out$grp %in% "A"]), 1, tolerance = 1e-8)
})

test_that("a group with zero mean fitness is named and gets NA relative fitness", {
  set.seed(3)
  df <- data.frame(
    surv = c(rbinom(30, 1, 0.6), rep(0, 20)),
    z = rnorm(50),
    yr = rep(c(2001, 2002), c(30, 20))
  )
  expect_warning(
    out <- suppressMessages(prepare_selection_data(df, "surv", "z", group = "yr")),
    "2002"
  )
  expect_true(all(is.na(out$relative_fitness[out$yr == 2002])))
  expect_false(anyNA(out$relative_fitness[out$yr == 2001]))
})

test_that("return_grouped skips rows whose group label is missing", {
  df <- sim_data(n = 120)
  df$grp[1:2] <- NA
  expect_warning(
    res <- suppressMessages(selection_coefficients(
      df, "w", c("z1", "z2"), fitness_type = "continuous", group = "grp", return_grouped = TRUE)),
    "missing"
  )
  expect_setequal(unique(res$Group), c("A", "B"))
})

test_that("a missing group label does not break the reference-group choice", {
  df <- sim_data(n = 120)
  df$year <- rep(c(2001, 2002), length.out = 120)
  df$year[120] <- NA
  prep <- suppressWarnings(suppressMessages(
    prepare_selection_data(df, "w", c("z1", "z2"), na_action = "none")))

  u <- suppressWarnings(suppressMessages(
    univariate_spline(prep, "w", "z1", fitness_type = "continuous", group = "year")))
  expect_s3_class(u, "univariate_fitness")
  expect_false(anyNA(u$grid$fit))

  s <- suppressWarnings(suppressMessages(
    correlated_fitness_surface(prep, "w", c("z1", "z2"), method = "gam", group = "year")))
  expect_false(anyNA(s$grid$.fit[s$grid$.inside]))
})

test_that("univariate_spline detects continuous fitness by default", {
  df <- sim_data(n = 80)
  df$z1 <- as.numeric(scale(df$z1))
  u <- suppressWarnings(suppressMessages(univariate_spline(df, "w", "z1")))
  expect_equal(u$fitness_type, "continuous")
  expect_match(u$ci_method, "parametric")
})

test_that("the p-value GLM has an intercept for each group", {
  set.seed(83)
  d <- data.frame(z1 = rnorm(160), z2 = rnorm(160), colony = rep(c("east", "west"), each = 80))
  d$alive <- rbinom(160, 1, plogis(ifelse(d$colony == "east", -1.2, 0.9) + 0.5 * d$z1))
  res <- quiet(selection_coefficients(d, "alive", c("z1", "z2"), fitness_type = "binary", group = "colony"))
  prep <- quiet(prepare_selection_data(d, "alive", c("z1", "z2"), group = "colony"))
  ref <- summary(glm(alive ~ colony + z1 + z2, family = binomial, data = prep))$coefficients
  lin <- res[res$Type == "Linear", ]
  expect_equal(lin$P_Value[match(c("z1", "z2"), lin$Term)], unname(ref[c("z1", "z2"), 4]))
})

test_that("a group with one level is fitted without a group term", {
  set.seed(41)
  d <- data.frame(z1 = rnorm(90), z2 = rnorm(90), site = "north")
  d$surv <- rbinom(90, 1, plogis(-0.5 + 0.8 * d$z1))
  prep <- quiet(prepare_selection_data(d, "surv", c("z1", "z2"), group = "site"))
  expect_message(sp <- suppressWarnings(univariate_spline(prep, "surv", "z1", group = "site")),
                 "one level")
  plain <- quiet(univariate_spline(prep, "surv", "z1"))
  expect_equal(fitted(sp$model), fitted(plain$model))
  expect_message(surf <- suppressWarnings(correlated_fitness_surface(prep, "surv", c("z1", "z2"), group = "site")),
                 "one level")
  expect_false(isTRUE(surf$group_effect))

  # a second site whose birds all lack the trait leaves one level in the fit
  south <- data.frame(z1 = NA_real_, z2 = rnorm(10), site = "south", surv = 1)
  both <- rbind(prep[, c("z1", "z2", "site", "surv")], south)
  expect_message(one <- suppressWarnings(univariate_spline(both, "surv", "z1", group = "site")), "one level")
  expect_equal(fitted(one$model), fitted(plain$model))
})

test_that("adaptive_landscape can use a surface fitted with a group term", {
  df <- sim_data(n = 120)
  prep <- suppressWarnings(suppressMessages(prepare_selection_data(df, "w", c("z1", "z2"), group = "grp")))
  surf <- suppressWarnings(suppressMessages(
    correlated_fitness_surface(prep, "w", c("z1", "z2"), method = "gam", group = "grp")))
  land <- suppressWarnings(suppressMessages(capture.output(
    al <- adaptive_landscape(prep, surf$model, c("z1", "z2"), grid_n = 6, simulation_n = 40))))
  expect_s3_class(al, "adaptive_landscape")
  expect_false(anyNA(al$grid$.mean_fit))
})

test_that("the 3D landscape is drawn from the grid without interpolation", {
  skip_if_not_installed("fields")
  g <- expand.grid(z1 = seq(-1, 1, length.out = 7), z2 = seq(-1, 1, length.out = 7))
  g$.mean_fit <- 1 - g$z1^2 - g$z2^2
  land <- structure(list(grid = g, trait_cols = c("z1", "z2")), class = "adaptive_landscape")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot_adaptive_landscape_3d(land, c("z1", "z2")))
})

test_that("the comparison overlay renders without optimum points", {
  g <- expand.grid(z1 = 1:5, z2 = 1:5)
  df <- rbind(
    cbind(g, fitness = runif(25), type = "Correlated Fitness (Individual)"),
    cbind(g, fitness = runif(25), type = "Adaptive Landscape (Population)")
  )
  plots <- plot_fitness_surfaces_comparison(list(combined_data = df, trait_cols = c("z1", "z2")))
  expect_no_error(ggplot2::ggplot_build(plots$overlay))
  plots2 <- plot_fitness_surfaces_comparison(
    list(combined_data = df, trait_cols = c("z1", "z2"),
         optimum_individual = data.frame(z1 = 2, z2 = 3),
         optimum_population = data.frame(z1 = 3, z2 = 2)),
    show_optima = FALSE)
  expect_no_error(ggplot2::ggplot_build(plots2$overlay))
})

test_that("selection_report differentials use the rows the gradients use", {
  df <- sim_data(n = 150)
  df$z2[1:30] <- NA # these individuals never reach the gradient models
  rep <- suppressWarnings(suppressMessages(
    selection_report(df, "w", c("z1", "z2"), fitness_type = "continuous")))

  cc <- df[stats::complete.cases(df[, c("w", "z1", "z2")]), ]
  z <- as.numeric(scale(cc$z1))
  w <- cc$w / mean(cc$w)
  S_ref <- mean((z - mean(z)) * (w - mean(w)))
  expect_equal(rep$Estimate[rep$Type == "Differential" & rep$Term == "z1"], S_ref, tolerance = 1e-10)
})

test_that("selection_report treats missing group labels as the gradients do", {
  df <- sim_data(n = 160, seed = 23)
  df$site <- rep(c("east", "west"), 80)
  df$site[c(3, 8, 20, 41, 77, 90, 101, 133)] <- NA
  labelled <- df
  labelled$site[is.na(labelled$site)] <- "none"
  a <- suppressWarnings(suppressMessages(selection_report(df, "w", c("z1", "z2"), fitness_type = "continuous", group = "site")))
  b <- suppressWarnings(suppressMessages(selection_report(labelled, "w", c("z1", "z2"), fitness_type = "continuous", group = "site")))
  # S and beta both count the unlabelled rows as one more group
  expect_equal(a$Estimate, b$Estimate)
})

test_that("selection_report labels the fitness type actually used", {
  d <- data.frame(surv = rbinom(80, 1, 0.5), z1 = rnorm(80), z2 = rnorm(80))
  rep <- suppressWarnings(suppressMessages(
    selection_report(d, "surv", c("z1", "z2"), fitness_type = "continuous")))
  expect_equal(attr(rep, "fitness_type"), "continuous")
  expect_output(print(rep), "Fitness type: continuous")
})

test_that("bootstrap_selection validates n_boot and resamples within groups", {
  df <- sim_data(n = 60)
  expect_error(bootstrap_selection(df, "w", c("z1", "z2"), n_boot = NA), "n_boot")
  expect_error(bootstrap_selection(df, "w", c("z1", "z2"), n_boot = 25.5), "n_boot")

  set.seed(2)
  b <- suppressWarnings(suppressMessages(
    bootstrap_selection(df, "w", c("z1", "z2"), fitness_type = "continuous",
                        group = "grp", n_boot = 30)))
  expect_true(all(b$N_Boot <= 30))
  ref <- suppressWarnings(suppressMessages(
    selection_coefficients(df, "w", c("z1", "z2"), fitness_type = "continuous", group = "grp")))
  expect_equal(b$Estimate, ref$Beta_Coefficient)
})

test_that("the assumption checks take the group out, as the gradient models do", {
  set.seed(83)
  d <- data.frame(z1 = rnorm(300), z2 = rnorm(300), yr = rep(c("y1", "y2", "y3"), each = 100))
  d$kids <- rpois(300, exp(c(y1 = -0.2, y2 = 1.2, y3 = 2.2)[d$yr] + 0.2 * d$z1 - 0.1 * d$z2^2))
  chk <- suppressWarnings(check_selection_assumptions(d, "kids", c("z1", "z2"), group = "yr"))
  # once each year has its own mean the counts are close to Poisson
  expect_lt(chk$statistic[grepl("dispersion", chk$check)], 1.5)
  rep <- suppressWarnings(suppressMessages(selection_report(d, "kids", c("z1", "z2"), group = "yr")))
  expect_equal(unname(attr(rep, "p_model")[1]), "poisson(log)")
})

test_that("a population where every bird lived is not taken for separation", {
  set.seed(9)
  d <- data.frame(z1 = rnorm(200), z2 = rnorm(200), pop = rep(c("p1", "p2"), c(170, 30)))
  d$alive <- rbinom(200, 1, plogis(0.2 + 0.7 * d$z1))
  d$alive[d$pop == "p2"] <- 1
  seen <- character()
  withCallingHandlers(suppressMessages(selection_coefficients(d, "alive", c("z1", "z2"), group = "pop")),
                      warning = function(w) {
                        seen <<- c(seen, conditionMessage(w))
                        invokeRestart("muffleWarning")
                      })
  expect_false(any(grepl("separation", seen)))
  # nor in the checks, where that population is fitted at 1 by its intercept
  chk <- suppressWarnings(check_selection_assumptions(d, "alive", c("z1", "z2"), group = "pop"))
  expect_equal(chk$note[grepl("separation", chk$check)], "no sign of separation")
})

test_that("the dispersion note names the model each set of gradients uses", {
  set.seed(5)
  d <- data.frame(z1 = rnorm(400), z2 = rnorm(400))
  d$kids <- rpois(400, exp(1.6 - 0.7 * d$z1^2 + 0.1 * d$z2))
  chk <- suppressWarnings(check_selection_assumptions(d, "kids", c("z1", "z2")))
  expect_match(chk$note[grepl("dispersion", chk$check)],
               "linear p-values use a negative binomial model, the quadratic ones a Poisson model")
  # unlabelled rows count as one more group
  set.seed(3)
  g <- data.frame(z1 = rnorm(360), z2 = rnorm(360), site = rep(c("a", "b", NA), each = 120))
  mu <- exp(1 + 0.2 * g$z1)
  g$kids <- ifelse(is.na(g$site), rnbinom(360, mu = mu, size = 0.8), rpois(360, mu))
  chk <- suppressWarnings(check_selection_assumptions(g, "kids", c("z1", "z2"), group = "site"))
  expect_equal(attr(chk, "n"), 360)
  rep <- suppressWarnings(suppressMessages(selection_report(g, "kids", c("z1", "z2"), group = "site")))
  expect_equal(unname(attr(rep, "p_model")), rep("negative binomial", 2))
  expect_match(chk$note[grepl("dispersion", chk$check)], "the p-values use a negative binomial model")
})

test_that("the spline and surface keep unlabelled rows as one more group and fit a numeric group by level", {
  set.seed(5)
  g <- data.frame(z1 = rnorm(300), z2 = rnorm(300), site = rep(c("a", "b", NA), each = 100))
  g$w <- 1 + 0.3 * g$z1 + ifelse(is.na(g$site), 0.5, 0) + rnorm(300, 0, 0.3)
  sp <- suppressWarnings(suppressMessages(univariate_spline(g, "w", "z1", fitness_type = "continuous", group = "site")))
  expect_equal(nrow(sp$model$model), 300)
  expect_true("siteNA" %in% names(coef(sp$model)))
  sf <- suppressWarnings(suppressMessages(correlated_fitness_surface(g, "w", c("z1", "z2"), method = "gam", group = "site")))
  expect_equal(nrow(sf$model$model), 300)
  expect_true("siteNA" %in% names(coef(sf$model)))
  g$year <- rep(c(2001, 2002, 2003), 100)
  sy <- suppressWarnings(suppressMessages(univariate_spline(g, "w", "z1", fitness_type = "continuous", group = "year")))
  expect_equal(grep("^year", names(coef(sy$model)), value = TRUE), c("year2002", "year2003"))
})

test_that("the canonical bootstrap resamples unlabelled rows as one more group", {
  set.seed(6)
  d <- data.frame(z1 = rnorm(120), z2 = rnorm(120), g = rep(c("a", "b", NA), each = 40))
  d$w <- 1 + 0.2 * d$z1 - 0.1 * d$z2^2 + rnorm(120, 0, 0.3)
  real <- selection_coefficients
  sizes <- integer()
  local_mocked_bindings(selection_coefficients = function(data, ...) {
    sizes <<- c(sizes, nrow(data))
    real(data, ...)
  })
  quiet(canonical_analysis(d, "w", c("z1", "z2"), group = "g", bootstrap = TRUE, n_boot = 5))
  expect_gt(length(sizes), 5)
  expect_true(all(sizes == 120))
})

test_that("Mardia's subsample is the same every run and leaves the random numbers alone", {
  set.seed(1)
  big <- data.frame(z1 = rnorm(2500), z2 = rexp(2500))
  big$w <- big$z1 + rnorm(2500)
  set.seed(42)
  ref <- runif(1)
  set.seed(42)
  a <- suppressWarnings(check_selection_assumptions(big, "w", c("z1", "z2")))
  expect_equal(runif(1), ref)
  set.seed(7)
  b <- suppressWarnings(check_selection_assumptions(big, "w", c("z1", "z2")))
  expect_equal(a$statistic, b$statistic)
  expect_match(a$note[1], "random 2000 of 2500 rows")
})

test_that("a numeric group is held at the same level on the surface and the landscape", {
  # level 1 is the most common, level 2 the middle one: the landscape used to
  # take the middle level of the raw column while the surface took the most
  # common level of the factor it fitted
  set.seed(11)
  n <- 300
  d <- data.frame(g = rep(c(1, 2, 3), c(180, 60, 60)), z1 = rnorm(n), z2 = rnorm(n))
  d$w <- rbinom(n, 1, plogis(-0.3 + 0.6 * d$z1 - 0.4 * d$z2^2 + 0.5 * (d$g == 3)))
  surf_msgs <- character()
  surf <- withCallingHandlers(
    correlated_fitness_surface(d, "w", c("z1", "z2"), group = "g", grid_n = 10),
    message = function(m) { surf_msgs <<- c(surf_msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  land_msgs <- capture_messages(adaptive_landscape(d, surf$model, c("z1", "z2"), simulation_n = 50, grid_n = 5))
  surf_ref <- sub(".*group = '([^']+)'.*", "\\1", grep("as reference", surf_msgs, value = TRUE))
  land_ref <- sub(".*reference level '([^']+)'.*", "\\1", grep("reference level", land_msgs, value = TRUE))
  expect_identical(surf_ref, "1")
  expect_identical(land_ref, surf_ref)
})
