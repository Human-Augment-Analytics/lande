# bootstrap_selection() refits only the least-squares gradients on each
# resample, with .ols_gradients(), when that gives the point estimate's
# gradients exactly, and otherwise all of selection_coefficients() as before.
# These tests hold the helper to the least-squares models fitted the long way
# (prepare_selection_data(), then the analyze_*() models), and the bootstrap
# to the loop it replaced, which ran selection_coefficients() on every
# resample. selection_coefficients() does not use the helper, so neither
# comparison checks the helper against itself.
#
# The comparisons are bit for bit, the proof that nothing moved. On CRAN they
# allow a relative difference of 1e-12 instead, so that a change in R's
# arithmetic that moves a last bit cannot fail the package's checks there.

# testthat's own rule for skip_on_cran(); test_local() and devtools set NOT_CRAN
on_cran <- function() {
  env <- Sys.getenv("NOT_CRAN")
  if (identical(env, "")) !interactive() else !isTRUE(as.logical(env))
}
expect_same_bits <- function(object, expected) {
  if (on_cran()) {
    expect_equal(object, expected, tolerance = 1e-12)
  } else {
    expect_identical(object, expected)
  }
}

# the gradients as the least-squares models give them: prepare, fit, extract
ls_reference <- function(data, fit, traits, fitness_type = "auto", standardize = TRUE,
                         group = NULL, use_relative_for_fit = TRUE) {
  rel <- paste0(fit, "_relative")
  prep <- prepare_selection_data(data, fit, traits, standardize = standardize, group = group,
                                 name_relative = rel)
  if (fitness_type == "auto") fitness_type <- detect_family(prep[[fit]])$type
  resp <- if (use_relative_for_fit) rel else fit
  raw <- if (fitness_type %in% c("binary", "count")) fit else NULL
  lin <- analyze_linear_selection(prep, resp, traits, fitness_type, binary_response_col = raw, group = group)
  nl <- analyze_nonlinear_selection(prep, resp, traits, fitness_type, binary_response_col = raw, group = group)
  rbind(extract_linear_coefficients(traits, lin),
        extract_quadratic_coefficients(traits, nl),
        extract_interaction_coefficients(traits, nl))
}

# the helper's gradients for the terms of a coefficient table
helper_gradients <- function(data, fit, traits, standardize = TRUE, group = NULL,
                             use_relative_for_fit = TRUE, terms) {
  cols <- lapply(stats::setNames(nm = unique(c(fit, traits))), function(v) data[[v]])
  g <- .ols_gradients(cols, fit, traits, .group_codes(data, group),
                      standardize = standardize, use_relative_for_fit = use_relative_for_fit)
  g[match(paste(terms$Term, terms$Type), .gradient_keys(traits))]
}

expect_helper_matches <- function(data, fit, traits, fitness_type = "auto", standardize = TRUE,
                                  group = NULL, use_relative_for_fit = TRUE) {
  ref <- quiet(ls_reference(data, fit, traits, fitness_type, standardize, group, use_relative_for_fit))
  expect_gt(nrow(ref), 0)
  expect_same_bits(helper_gradients(data, fit, traits, standardize, group, use_relative_for_fit, ref),
                   ref$Beta_Coefficient)
}

test_that("the gradient helper gives the least-squares models' gradients bit for bit", {
  tr <- c("weight", "total_length")
  expect_helper_matches(bumpus, "survival", tr, "binary", group = "sex")
  expect_helper_matches(bumpus, "survival", tr, "binary")
  expect_helper_matches(bumpus, "survival", c(tr, "humerus"), "binary", group = "sex")
  expect_helper_matches(bumpus, "survival", "weight", "binary")
  expect_helper_matches(finch_community, "lifespan", c("beak_length", "beak_depth"), "count")

  cp <- crescent_pond_pupfish
  cp <- cp[stats::complete.cases(cp[, c("ln.growth", "jaw", "body", "eye")]), ]
  cp$growth <- exp(cp$ln.growth)
  expect_helper_matches(cp, "growth", c("jaw", "body", "eye"), "continuous", group = "density")

  # integer traits, unstandardised: their products must not overflow
  bi <- bumpus
  bi$w_int <- as.integer(round(bi$weight * 1e4))
  bi$l_int <- as.integer(round(bi$total_length * 1e3))
  expect_true(is.na(suppressWarnings(bi$w_int[1] * bi$l_int[1])))
  expect_helper_matches(bi, "survival", c("w_int", "l_int"), standardize = FALSE)
  # and standardised, with a missing value (integers averaged as doubles then)
  bi$w_int[5] <- NA
  expect_helper_matches(bi, "survival", c("w_int", "l_int"), group = "sex")
  # with none missing, prepare_selection_data() averages an integer trait as
  # integers, which can differ in the last bit from averaging it as doubles;
  # on 64-bit ARM it does for this trait, pooled and in group "a"
  set.seed(4)
  d <- data.frame(x = as.integer(round(rnorm(200, 2e5, 5e4))), g = rep(c("a", "b"), each = 100))
  d$y <- rnorm(200)
  d$w <- rbinom(200, 1, 0.5)
  expect_helper_matches(d, "w", c("x", "y"), "binary")
  expect_helper_matches(d, "w", c("x", "y"), "binary", group = "g")
  # and with one missing it averages them as doubles, which differs here too
  d$x[3] <- NA
  expect_helper_matches(d, "w", c("x", "y"), "binary")

  # absolute fitness
  expect_helper_matches(bumpus, "survival", tr, "binary", group = "sex", use_relative_for_fit = FALSE)
})

test_that("the gradient helper copies the preparation's odd cases", {
  tr <- c("weight", "total_length")
  b <- bumpus
  b$weight[c(3, 10, 50)] <- NA
  b$sex[c(2, 30, 80)] <- NA
  expect_helper_matches(b, "survival", tr, "binary", group = "sex")

  # a single-bird group: its traits are only centred, and the helper says so
  b <- bumpus
  b$grp <- as.character(b$sex)
  b$grp[5] <- "solo"
  expect_helper_matches(b, "survival", tr, "binary", group = "grp")
  cols <- list(survival = b$survival, weight = b$weight, total_length = b$total_length)
  g <- .ols_gradients(cols, "survival", tr, .group_codes(b, "grp"))
  expect_true(attr(g, "degenerate"))
  expect_null(.ols_gradients(cols, "survival", tr, .group_codes(b, "grp"), drop_degenerate = TRUE))

  # a group where every bird died has no relative fitness and drops out
  b <- bumpus
  b$grp <- as.character(b$sex)
  b$grp[which(b$survival == 0)[1:6]] <- "dead"
  expect_helper_matches(b, "survival", tr, "binary", group = "grp")

  # an aliased trait: its terms are NA in the helper and absent from the table
  b <- bumpus
  b$w2 <- 2 * b$weight
  expect_helper_matches(b, "survival", c("weight", "w2"), "binary")
})

test_that("the gradient helper matches on resampled data, as the bootstrap uses it", {
  b <- bumpus
  b$grp <- as.character(b$sex)
  b$grp[c(1, 2, 60)] <- "tiny"
  tr <- c("weight", "total_length")
  set.seed(5)
  n_degenerate <- 0
  for (i in 1:25) {
    d <- b[sample.int(nrow(b), replace = TRUE), ]
    cols <- lapply(stats::setNames(nm = c("survival", tr)), function(v) d[[v]])
    g <- .ols_gradients(cols, "survival", tr, .group_codes(d, "grp"))
    n_degenerate <- n_degenerate + attr(g, "degenerate")
    ref <- quiet(ls_reference(d, "survival", tr, "binary", group = "grp"))
    expect_same_bits(g[match(paste(ref$Term, ref$Type), .gradient_keys(tr))], ref$Beta_Coefficient)
  }
  # some of those resamples lost a trait's variance in the tiny group
  expect_gt(n_degenerate, 0)
})

test_that("the gradient helper names the terms as selection_coefficients() does", {
  tr <- c("weight", "total_length", "humerus")
  sc <- quiet(selection_coefficients(bumpus, "survival", tr, fitness_type = "binary"))
  expect_identical(.gradient_keys(tr), paste(sc$Term, sc$Type))
  expect_identical(.gradient_keys("weight"), c("weight Linear", "weight\u00b2 Quadratic"))
})

# The bootstrap as it was: selection_coefficients() on every resample, kept
# only for its gradients.
reference_bootstrap <- function(data, fitness_col, trait_cols, fitness_type = "auto",
                                standardize = TRUE, group = NULL, use_relative_for_fit = TRUE,
                                n_boot = 1000, conf = 0.95) {
  fit_once <- function(d) {
    selection_coefficients(d, fitness_col, trait_cols, fitness_type = fitness_type,
                           standardize = standardize, group = group,
                           use_relative_for_fit = use_relative_for_fit)
  }
  point <- suppressMessages(fit_once(data))
  keys <- paste(point$Term, point$Type)
  fit_resample <- function(d) {
    degenerate <- FALSE
    res <- tryCatch(
      withCallingHandlers(
        suppressMessages(fit_once(d)),
        warning = function(w) {
          if (grepl("no variance|Zero-variance", conditionMessage(w))) degenerate <<- TRUE
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) NULL
    )
    if (degenerate) NULL else res
  }
  n <- nrow(data)
  strata <- if (is.null(group)) rep(1L, n) else as.integer(addNA(factor(data[[group]]), ifany = TRUE))
  strata_rows <- split(seq_len(n), strata, drop = TRUE)
  resample_rows <- function() {
    unlist(lapply(strata_rows, function(r) r[sample.int(length(r), length(r), replace = TRUE)]),
           use.names = FALSE)
  }
  boot <- matrix(NA_real_, nrow = length(keys), ncol = n_boot, dimnames = list(keys, NULL))
  for (b in seq_len(n_boot)) {
    res <- fit_resample(data[resample_rows(), , drop = FALSE])
    if (is.null(res)) next
    boot[, b] <- res$Beta_Coefficient[match(keys, paste(res$Term, res$Type))]
  }
  n_used <- rowSums(!is.na(boot))
  n_min <- min(n_used)
  if (n_min < 2) stop("Bootstrap failed: fewer than 2 usable resamples for at least one term")
  if (n_min < n_boot) {
    warning(n_boot - n_min, " of ", n_boot,
            " resamples were unusable for at least one term and were dropped")
  }
  a <- (1 - conf) / 2
  out <- data.frame(
    Term = point$Term, Type = point$Type, Estimate = point$Beta_Coefficient,
    Boot_SE = apply(boot, 1, stats::sd, na.rm = TRUE),
    CI_lower = apply(boot, 1, stats::quantile, probs = a, na.rm = TRUE),
    CI_upper = apply(boot, 1, stats::quantile, probs = 1 - a, na.rm = TRUE),
    P_Value = point$P_Value, N_Boot = as.integer(n_used),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  rownames(out) <- NULL
  attr(out, "n_boot") <- n_min
  attr(out, "conf") <- conf
  out
}

# run a bootstrap from a fixed seed: its result, its warnings and the RNG state after
seeded_run <- function(FUN, ..., seed = 42) {
  warns <- character()
  set.seed(seed)
  value <- withCallingHandlers(
    suppressMessages(FUN(...)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warns, seed = get(".Random.seed", envir = globalenv()))
}

expect_same_bootstrap <- function(...) {
  new <- seeded_run(bootstrap_selection, ...)
  old <- seeded_run(reference_bootstrap, ...)
  expect_same_bits(new$value, old$value)
  expect_identical(new$warnings, old$warnings)
  expect_identical(new$seed, old$seed)
  invisible(new)
}

# counts the calls bootstrap_selection() makes to selection_coefficients() and
# to the p-value GLM while running `code`
count_fits <- function(code) {
  n <- c(selection_coefficients = 0L, glm = 0L)
  sc <- selection_coefficients
  glm_fit <- .fit_pvalue_glm
  local_mocked_bindings(
    selection_coefficients = function(...) {
      n[["selection_coefficients"]] <<- n[["selection_coefficients"]] + 1L
      sc(...)
    },
    .fit_pvalue_glm = function(...) {
      n[["glm"]] <<- n[["glm"]] + 1L
      glm_fit(...)
    }
  )
  quiet(code)
  n
}

test_that("the bootstrap equals selection_coefficients() refitted on the same resamples", {
  tr <- c("weight", "total_length")
  expect_same_bootstrap(bumpus, "survival", tr, fitness_type = "binary", group = "sex", n_boot = 30)
  expect_same_bootstrap(bumpus, "survival", c(tr, "humerus"), fitness_type = "binary", n_boot = 20)
  expect_same_bootstrap(finch_community, "lifespan", c("beak_length", "beak_depth"),
                        fitness_type = "count", n_boot = 6)

  cp <- crescent_pond_pupfish
  cp <- cp[stats::complete.cases(cp[, c("ln.growth", "jaw", "body")]), ]
  cp$growth <- exp(cp$ln.growth)
  expect_same_bootstrap(cp, "growth", c("jaw", "body"), fitness_type = "continuous",
                        group = "density", n_boot = 12, conf = 0.9)

  bi <- bumpus
  bi$w_int <- as.integer(round(bi$weight * 1e4))
  bi$l_int <- as.integer(round(bi$total_length * 1e3))
  expect_same_bootstrap(bi, "survival", c("w_int", "l_int"), standardize = FALSE, n_boot = 20)
  expect_same_bootstrap(bumpus, "survival", tr, group = "sex", use_relative_for_fit = FALSE, n_boot = 20)

  b <- bumpus
  b$weight[c(3, 10, 50)] <- NA
  b$sex[c(2, 30, 80)] <- NA
  expect_same_bootstrap(b, "survival", tr, fitness_type = "binary", group = "sex", n_boot = 30)
})

test_that("the bootstrap drops the same degenerate resamples, with the same warning", {
  b <- bumpus
  b$grp <- as.character(b$sex)
  b$grp[c(which(b$survival == 1)[1], which(b$survival == 0)[1:2])] <- "tiny"
  run <- expect_same_bootstrap(b, "survival", c("weight", "total_length"),
                               fitness_type = "binary", group = "grp", n_boot = 60)
  expect_lt(attr(run$value, "n_boot"), 60)
  expect_match(run$warnings, "resamples were unusable", all = FALSE)
})

test_that("resamples whose least-squares fit is exact are kept, as the old loop now keeps them", {
  # 4 complete birds of 8: a resample with 3 distinct complete birds fits the
  # quadratic model (3 coefficients) exactly
  set.seed(1)
  tiny <- data.frame(w = rexp(8), x = rnorm(8))
  tiny$x[1:4] <- NA
  run <- expect_same_bootstrap(tiny, "w", "x", fitness_type = "continuous", n_boot = 40)
  # before .compute_vif() skipped the NaN VIFs of an exact fit, its check
  # stopped selection_coefficients() on them, and the old loop dropped them
  skip_if_not_installed("car")
  vif_was <- function(fit) {
    if (length(attr(stats::terms(fit), "term.labels")) < 2) return(NULL)
    v <- tryCatch(car::vif(fit), error = function(e) {
      warning("VIF calculation failed: ", e$message)
      NULL
    })
    if (!is.null(v) && any(v > 5)) warning("Collinear traits (VIF above 5) may inflate the standard errors")
    v
  }
  local_mocked_bindings(.compute_vif = vif_was)
  before <- seeded_run(reference_bootstrap, tiny, "w", "x", fitness_type = "continuous", n_boot = 40)
  expect_lt(attr(before$value, "n_boot"), attr(run$value, "n_boot"))
})

test_that("the bootstrap refits all of selection_coefficients() when the shortcut misses the point estimate", {
  # the preparation overwrites a trait named like the relative-fitness column,
  # and the model formula reads a trait named "log(weight)" as log() of weight;
  # the shortcut fits the columns themselves, so it falls back to the old loop
  b <- bumpus
  b$survival_relative <- b$weight
  b[["log(weight)"]] <- log(b$weight)
  for (tr in list(c("survival_relative", "total_length"), c("log(weight)", "total_length"))) {
    expect_same_bootstrap(b, "survival", tr, fitness_type = "binary", n_boot = 20)
    n <- count_fits(bootstrap_selection(b, "survival", tr, fitness_type = "binary", n_boot = 20))
    expect_identical(n[["selection_coefficients"]], 21L)
  }
  # the old loop's drop rule comes with it: a tiny group's degenerate resamples
  b$grp <- as.character(b$sex)
  b$grp[c(which(b$survival == 1)[1], which(b$survival == 0)[1:2])] <- "tiny"
  run <- expect_same_bootstrap(b, "survival", c("log(weight)", "total_length"),
                               fitness_type = "binary", group = "grp", n_boot = 40)
  expect_lt(attr(run$value, "n_boot"), 40)
})

test_that("ordinary data take the shortcut, and fit the p-value GLMs only on the full data", {
  # the shortcut needs the helper to give lm()'s gradients bit for bit
  skip_on_cran()
  n <- count_fits(bootstrap_selection(bumpus, "survival", c("weight", "total_length"),
                                      group = "sex", n_boot = 20))
  # one point estimate, whose linear and quadratic models each have a GLM
  expect_identical(n, c(selection_coefficients = 1L, glm = 2L))
  n <- count_fits(bootstrap_selection(finch_community, "lifespan", c("beak_length", "beak_depth"),
                                      n_boot = 5))
  expect_identical(n, c(selection_coefficients = 1L, glm = 2L))
  b <- bumpus
  b$weight[c(3, 10, 50)] <- NA
  b$sex[c(2, 30, 80)] <- NA
  n <- count_fits(bootstrap_selection(b, "survival", c("weight", "total_length", "humerus"),
                                      group = "sex", n_boot = 20))
  expect_identical(n[["selection_coefficients"]], 1L)
})

test_that("a resample whose p-value GLM fails is kept, where the old loop dropped it", {
  # the GLM's failure is numerical
  skip_on_cran()
  fc <- finch_community
  fc <- fc[stats::complete.cases(fc[, c("lifespan", "beak_length", "beak_depth")]), ][1:120, ]
  tr <- c("beak_length", "beak_depth")
  new <- seeded_run(bootstrap_selection, fc, "lifespan", tr, n_boot = 15, seed = 1)
  old <- seeded_run(reference_bootstrap, fc, "lifespan", tr, n_boot = 15, seed = 1)
  expect_equal(attr(new$value, "n_boot"), 15)
  expect_true(all(new$value$N_Boot == 15L))
  expect_false(any(grepl("unusable", new$warnings)))
  # one of the old loop's resamples failed in the quadratic count GLM
  expect_equal(attr(old$value, "n_boot"), 14)
  same <- c("Term", "Type", "Estimate", "P_Value")
  expect_identical(new$value[same], old$value[same])
  expect_identical(new$seed, old$seed)
})

test_that("the data's own <fitness>_relative column does not reach the resamples", {
  # 30 birds, one survivor, and an unrelated column named survival_relative
  set.seed(123)
  sm <- bumpus[c(which(bumpus$survival == 0)[1:29], which(bumpus$survival == 1)[1]),
               c("survival", "weight", "total_length")]
  sm$survival_relative <- round(runif(nrow(sm)), 3)
  tr <- c("weight", "total_length")
  with_col <- seeded_run(bootstrap_selection, sm, "survival", tr, n_boot = 60, seed = 61)
  without <- seeded_run(bootstrap_selection, sm[c("survival", tr)], "survival", tr, n_boot = 60, seed = 61)
  expect_identical(with_col, without)
  # resamples without the survivor have no relative fitness and are dropped
  expect_lt(attr(with_col$value, "n_boot"), 60)
  # the old loop fitted the data's own column on them instead
  skip_on_cran()
  old <- seeded_run(reference_bootstrap, sm, "survival", tr, n_boot = 60, seed = 61)
  expect_equal(attr(old$value, "n_boot"), 60)
})
