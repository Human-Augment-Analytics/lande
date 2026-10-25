# ======================================================
# univariate_spline.R
# Estimate univariate correlated fitness function
#
# Important Concept Explanation:
# This function calculates a UNIVARIATE CORRELATED FITNESS FUNCTION
# Definition: Individual fitness ~ Individual phenotype (single trait)
# Formula: w ~ f(z)
#
# This is a special case of correlated fitness surface (1D instead of 2D)
# Adaptive landscape would require: Mean fitness ~ Population mean phenotype
#
# IMPORTANT NOTES
#   - Traits should be standardized BEFORE calling this function
#   - Use prepare_selection_data() with standardize = TRUE and optional group
#   - The GAM smooth term estimates the shape of the fitness function
#   - DO NOT use scale() within this function (double standardization)
# ======================================================


#' @noRd
# internal utility: one independent fit per level of a grouping column. A group
# whose fit fails is left out with a warning that names it.
.fit_by_group <- function(data, group, fit) {
  if (is.null(group)) stop("`by_group = TRUE` needs a `group` column")
  if (!group %in% names(data)) stop("Group column '", group, "' not found in data")
  g <- data[[group]]
  seen <- g[!is.na(g)]
  levels_used <- if (is.factor(seen)) levels(droplevels(seen)) else sort(unique(seen))
  fits <- lapply(levels_used, function(level) {
    rows <- data[!is.na(g) & g == level, , drop = FALSE]
    tryCatch(fit(rows), error = function(e) {
      warning("Group '", level, "' left out: ", conditionMessage(e), call. = FALSE)
      NULL
    })
  })
  names(fits) <- as.character(levels_used)
  fits <- fits[!vapply(fits, is.null, logical(1))]
  if (!length(fits)) stop("No group could be fitted")
  attr(fits, "group") <- group
  fits
}

#' @noRd
# internal utility: the options a refit on a worker is given from the
# caller's session. The contrasts and matprod can change the fit; the others
# change which warnings a fit gives or how they are worded.
.refit_options <- function() {
  nms <- c("contrasts", "matprod", "digits", "OutDec", "scipen",
           "warnPartialMatchArgs", "warnPartialMatchAttr", "warnPartialMatchDollar")
  stats::setNames(lapply(nms, getOption), nms)
}

#' @noRd
# internal utility: one bootstrap refit of the spline on rows i of d, run on a
# worker of a parallel cluster. Returns the curve over the grid on the
# response scale, whether the fit failed, the warnings and messages of the
# fit and of the prediction apart, and an error in the prediction, which the
# serial loop does not catch; the caller raises them where that loop would.
# The worker fits under the caller's random number generator kinds, as mgcv
# draws from it (and puts it back) when a thin-plate basis has to choose its
# knots, and under the caller's options from .refit_options(). Only base R,
# stats and mgcv are called, so it runs on a worker without lande.
.spline_refit <- function(i, d, fml, fam, smoothing, grid, rng_kind, opts) {
  kind <- RNGkind()
  if (!identical(kind, rng_kind)) {
    try(suppressWarnings(RNGkind(rng_kind[1], rng_kind[2], rng_kind[3])), silent = TRUE)
    on.exit(try(suppressWarnings(RNGkind(kind[1], kind[2], kind[3])), silent = TRUE), add = TRUE)
  }
  old <- options(opts)
  on.exit(options(old), add = TRUE)
  said <- list()
  keep <- function(cnd) {
    said[[length(said) + 1L]] <<- cnd
    tryInvokeRestart(if (inherits(cnd, "warning")) "muffleWarning" else "muffleMessage")
  }
  fit_b <- withCallingHandlers(
    tryCatch(
      mgcv::gam(fml, data = d[i, , drop = FALSE], family = fam, method = smoothing, na.action = stats::na.omit),
      error = function(e) NULL
    ),
    warning = keep, message = keep
  )
  out <- list(fit_said = said, failed = is.null(fit_b))
  if (out$failed) return(out)
  said <- list()
  err <- NULL
  out$curve <- tryCatch(
    withCallingHandlers(
      fit_b$family$linkinv(as.numeric(stats::predict(fit_b, newdata = grid, type = "link"))),
      warning = keep, message = keep
    ),
    error = function(e) {
      err <<- e
      NULL
    }
  )
  out$predict_said <- said
  out$error <- err
  out
}

#' Estimate univariate correlated fitness function
#'
#' This function calculates a univariate correlated fitness function using a GAM smooth term.
#' The formula is \code{w ~ f(z)}. Traits should be standardized BEFORE calling this function.
#'
#' @param data A data frame containing fitness and trait measurements.
#' @param fitness_col A string specifying the name of the fitness column.
#' @param trait_col A string specifying the name of the trait column (must be numeric and standardized).
#' @param fitness_type A string indicating the fitness type: \code{"auto"} (detect from the data, the default), \code{"binary"}, \code{"count"}, or \code{"continuous"}. Binary and count fitness are fitted on the raw values with a binomial or Poisson family; continuous fitness on relative fitness with a Gaussian one.
#' @param group Optional string specifying a grouping variable. If provided, group fixed effects are included (none when it has only one level); rows with no group label are one more level.
#' @param relative_col Optional string naming a pre-computed relative fitness column to use for continuous fitness (e.g. one produced within groups by \code{prepare_selection_data}).
#' @param k Integer specifying the basis dimension for the smooth term. Default is 10.
#' @param bs Spline basis: \code{"cr"} (cubic regression spline, the default),
#'   \code{"tp"} (thin plate) or \code{"ps"} (P-spline).
#' @param smoothing How the smoothing parameter is chosen: \code{"GCV.Cp"}
#'   (generalised cross-validation, the default, which mgcv applies in its UBRE
#'   form to binary and Poisson fitness and as plain GCV to quasi-Poisson; a
#'   negative binomial fit uses REML in its place), \code{"REML"} or
#'   \code{"ML"}.
#'   The result's \code{spline_type} names the criterion mgcv used.
#' @param bootstrap Logical; if \code{TRUE} the 95\% ribbon is obtained by resampling individuals and refitting (Schluter 1988). The default \code{FALSE} uses the parametric Wald interval, which is fast; the bootstrap refits the spline \code{n_boot} times.
#' @param n_boot Integer number of bootstrap resamples used when \code{bootstrap = TRUE}. Default is 1000.
#' @param by_group Logical; if \code{TRUE} fit each level of \code{group} on its
#'   own rows and return a named list of fits, one per group.
#'   The default \code{FALSE} fits one curve with the group as a fixed
#'   effect. Prepare the data with the same \code{group} first, so that each
#'   group is standardised on its own.
#' @param count_family Family for count fitness: \code{"poisson"} (the
#'   default), \code{"quasipoisson"} or \code{"nb"}, as in
#'   \code{correlated_fitness_surface()}; the smoothing parameter is chosen
#'   under whichever is used. The result's \code{dispersion} is the Pearson
#'   dispersion of the fit, and a Poisson fit warns when it is above 1.5.
#' @param cl Optional cluster from \code{parallel::makeCluster()} to run the
#'   bootstrap refits on. The default \code{NULL} runs them one after another
#'   in this session. With a cluster the resamples are still drawn here, in
#'   the same order, and only the refits go to the workers, which fit under
#'   this session's random number generator kinds and its options for
#'   contrasts, matrix products and warnings. Their warnings and messages are
#'   raised here in the order the serial loop raises them, so the ribbon, the
#'   warnings, \code{options(warn = 2)} and the state of the random number
#'   generator afterwards (or when a handler stops the call at a warning) are
#'   the same as without a cluster. A negative binomial fit always runs here,
#'   as each refit starts from the theta of the one before. The caller starts
#'   and stops the cluster. The workers fit with mgcv, so they need the same
#'   version of it as this session on their library paths (for example
#'   \code{parallel::clusterCall(cl, function(p) .libPaths(p), .libPaths())});
#'   lande itself need not be loaded on them.
#'
#' @details By default the fitness function is a penalised cubic regression
#'   spline with the smoothing parameter chosen by generalised
#'   cross-validation (UBRE for binary and Poisson fitness), following Schluter
#'   (1988). \code{bs} and \code{smoothing}
#'   set the family of curves and the criterion, to match the smoother of
#'   another study; the amount of smoothing is still estimated from the data,
#'   though the criterion can change it a good deal.
#'   If mgcv's check suggests the basis dimension was too small a warning
#'   says so; the check works on the residuals, which for survival carry too
#'   little to test, so it is skipped there. With \code{bootstrap = TRUE} the
#'   result depends on the random seed; call \code{set.seed()} first for a
#'   reproducible ribbon.
#'
#' @return A list of class \code{"univariate_fitness"} containing the fitted GAM model, a prediction grid, and metadata.
#' @export
#'
#' @examples
#' prep <- prepare_selection_data(bumpus, "survival", "total_length")
#' uni <- univariate_spline(prep, "survival", "total_length")
#' plot_univariate_fitness(uni, "total_length")
#'
#' # a thin-plate basis with REML, to match another study's smoother
#' univariate_spline(prep, "survival", "total_length", bs = "tp", smoothing = "REML")$spline_type
univariate_spline <- function(data,
                              fitness_col,
                              trait_col,
                              fitness_type = c("auto", "binary", "count", "continuous"),
                              group = NULL,
                              relative_col = NULL,
                              k = 10,
                              bs = c("cr", "tp", "ps"),
                              smoothing = c("GCV.Cp", "REML", "ML"),
                              bootstrap = FALSE,
                              n_boot = 1000,
                              by_group = FALSE,
                              count_family = c("poisson", "quasipoisson", "nb"),
                              cl = NULL) {
  fitness_type <- match.arg(fitness_type)
  bs <- match.arg(bs)
  smoothing <- match.arg(smoothing)
  count_family <- match.arg(count_family)
  if (!is.null(cl) && !inherits(cl, "cluster")) {
    stop("`cl` must be a cluster from parallel::makeCluster() or NULL")
  }

  if (isTRUE(by_group)) {
    fits <- .fit_by_group(data, group, function(rows) {
      univariate_spline(rows, fitness_col, trait_col, fitness_type = fitness_type, group = NULL,
                        relative_col = relative_col, k = k, bs = bs, smoothing = smoothing,
                        bootstrap = bootstrap, n_boot = n_boot, count_family = count_family, cl = cl)
    })
    return(fits)
  }

  # Input validation
  if (length(trait_col) != 1L) {
    stop("`trait_col` must be a single column name.")
  }
  if (!trait_col %in% names(data)) {
    stop("`trait_col` not found in `data`.")
  }
  if (!is.numeric(data[[trait_col]])) {
    stop("`trait_col` must be numeric")
  }
  if (!fitness_col %in% names(data)) {
    stop("Fitness column '", fitness_col, "' not found in data")
  }
  if (fitness_type == "auto") {
    detected <- detect_family(data[[fitness_col]])$type
    fitness_type <- if (detected %in% c("binary", "count")) detected else "continuous"
    message("Fitness type detected: ", fitness_type)
  }

  # Check if group column exists
  if (!is.null(group) && !group %in% names(data)) {
    stop("Group column '", group, "' not found in data")
  }

  # the trait comes in standardised; warn when it does not look it
  z_mean <- mean(data[[trait_col]], na.rm = TRUE)
  z_sd <- sd(data[[trait_col]], na.rm = TRUE)

  if (abs(z_mean) > 0.1 || abs(z_sd - 1) > 0.1) {
    warning(
      "Trait '", trait_col, "' does not look standardised ",
      "(mean ", round(z_mean, 3), ", SD ", round(z_sd, 3), "), so it is used in its own units; ",
      "prepare_selection_data() standardises it"
    )
  }

  if (fitness_type == "continuous") {
    # For continuous fitness, use relative fitness. Prefer an explicitly named
    # column (e.g. one relativised within groups by prepare_selection_data).
    rel_source <- relative_col
    if (is.null(rel_source) && "relative_fitness" %in% names(data)) {
      rel_source <- "relative_fitness"
    }

    if (!is.null(rel_source)) {
      if (!rel_source %in% names(data)) {
        stop("relative_col '", rel_source, "' not found in data")
      }
      y <- data[[rel_source]]
      fit_note <- paste0("Using relative fitness column '", rel_source, "'")
    } else {
      # Compute relative fitness on the fly (warning: not group-specific)
      if (!is.null(group)) {
        warning(
          "A group is set but there is no relative fitness column; ",
          "run prepare_selection_data() first or pass relative_col."
        )
      }
      y <- data[[fitness_col]] / mean(data[[fitness_col]], na.rm = TRUE)
      fit_note <- "Computed relative fitness on the fly (pooled)"
    }
    fam <- stats::gaussian()
    family_name <- "gaussian"
  } else if (fitness_type == "count") {
    # Count fitness: Poisson (or count_family) on the raw counts, plotted on
    # the response scale
    y <- data[[fitness_col]]
    fam <- .count_family(count_family)
    family_name <- switch(count_family, poisson = "poisson(log)", quasipoisson = "quasipoisson(log)",
                          nb = "negative binomial")
    if (!.is_raw_fitness(y[!is.na(y)], "count")) {
      warning(
        "fitness_type = 'count' but the values are not all non-negative integers; ",
        "fitting them as counts anyway"
      )
    }
    fit_note <- "Using original counts"
  } else {
    # Binary fitness - use original 0/1 values
    y <- data[[fitness_col]]
    fam <- stats::binomial("logit")
    family_name <- "binomial(logit)"

    # Check if really binary
    unique_vals <- unique(y[!is.na(y)])
    if (!all(unique_vals %in% c(0, 1))) {
      warning("fitness_type = 'binary' but the values are not all 0/1; fitting them as binary anyway")
    }
    fit_note <- "Using original binary fitness"
  }

  df <- data
  df[[".y"]] <- y

  # Adjust k based on unique values (avoid GAM errors)
  n_unique <- length(unique(df[[trait_col]][complete.cases(df[[trait_col]])]))
  n_obs <- sum(complete.cases(df[[trait_col]], y))

  # k should not exceed n_unique - 1
  k_adj <- min(k, max(3, n_unique - 1))
  if (k_adj < k) {
    warning(
      "Reducing k from ", k, " to ", k_adj,
      " (only ", n_unique, " unique values)"
    )
    k <- k_adj
  }

  # Check if sample size is sufficient
  if (n_obs < k * 2) {
    warning(
      "Sample size (", n_obs, ") may be insufficient for k = ", k,
      ". Consider reducing k."
    )
  }

  # Default: cubic regression spline with GCV smoothing (Schluter 1988). bs = "cr"
  # is a cubic regression spline basis; mgcv's default is thin plate.
  smooth <- paste0("s(", trait_col, ", bs = '", bs, "', k = ", k, ")")
  # the group enters as a factor; rows with no label are one more level
  if (!is.null(group)) df[[group]] <- droplevels(addNA(factor(df[[group]]), ifany = TRUE))
  # levels are counted on the rows the model will use
  if (!is.null(group) && length(unique(df[[group]][complete.cases(df[[trait_col]], y, df[[group]])])) < 2) {
    message("No group term, as '", group, "' has only one level")
    group <- NULL
  }
  if (!is.null(group)) {
    fml <- stats::as.formula(paste0(".y ~ ", group, " + ", smooth))
    message("Including group fixed effect: '", group, "'")
  } else {
    fml <- stats::as.formula(paste0(".y ~ ", smooth))
  }

  last_error <- NULL
  fit_gam <- function(d) {
    tryCatch(
      mgcv::gam(fml, data = d, family = fam, method = smoothing, na.action = stats::na.omit),
      error = function(e) {
        last_error <<- conditionMessage(e)
        NULL
      }
    )
  }

  fit <- fit_gam(df)
  if (is.null(fit)) {
    stop("GAM fitting failed: ", last_error, "\nTry reducing k (currently k = ", k, ").")
  }
  if (!is.null(fit$converged) && !fit$converged) {
    warning("GAM algorithm did not fully converge")
  }
  dispersion <- if (fitness_type == "count") .check_dispersion(fit, count_family) else NULL

  # mgcv's test of whether the basis had room to bend: a low k-index with a
  # small p-value means the curve may look straighter than the data are. The
  # test works on the residuals, which for 0/1 fitness carry very little to
  # check (as mgcv's gam.check help says), so it is skipped for survival
  kc <- if (fitness_type != "binary") tryCatch(mgcv::k.check(fit), error = function(e) NULL)
  if (!is.null(kc) && any(kc[, "k-index"] < 1 & kc[, "p-value"] < 0.05, na.rm = TRUE)) {
    warning("k = ", k, " may be too small for '", trait_col, "' (mgcv k-index ",
            round(min(kc[, "k-index"]), 2), "); try a larger k")
  }

  # Create prediction grid across observed trait range
  rng <- range(df[[trait_col]], na.rm = TRUE)
  grid <- data.frame(seq(rng[1], rng[2], length.out = 200))
  names(grid) <- trait_col

  # If group was specified, predictions need a reference group
  if (!is.null(group)) {
    ref_group <- .reference_group(df[[group]])
    grid[[group]] <- ref_group
    message("Predictions use group = '", ref_group, "' as reference")
  }

  linkinv <- fit$family$linkinv
  grid$fit <- linkinv(as.numeric(stats::predict(fit, newdata = grid, type = "link")))

  # 95% ribbon: bootstrap individuals (Schluter 1988) or parametric Wald interval
  ci_method <- if (bootstrap) "bootstrap (percentile)" else "parametric (Wald)"
  if (bootstrap) {
    n <- nrow(df)
    boot_fits <- matrix(NA_real_, nrow = nrow(grid), ncol = n_boot)
    # mgcv's negative binomial family keeps theta between fits, so its refits
    # depend on their order and stay in this session
    if (!is.null(cl) && !(fitness_type == "count" && count_family == "nb")) {
      # the resamples drawn as the loop below draws them, keeping the random
      # number state after each, then refitted on the cluster. The refit and
      # the formula travel without their environments (lande's namespace, and
      # this call's, which holds the data and every resample), and only the
      # model's columns go with them
      env <- globalenv()
      idx <- seeds <- vector("list", n_boot)
      for (b in seq_len(n_boot)) {
        idx[[b]] <- sample.int(n, n, replace = TRUE)
        seeds[b] <- list(get0(".Random.seed", envir = env, inherits = FALSE))
      }
      refit <- .spline_refit
      environment(refit) <- baseenv()
      fml_w <- fml
      environment(fml_w) <- baseenv()
      refits <- parallel::parLapply(cl, idx, refit,
                                    d = df[, unique(c(".y", trait_col, group)), drop = FALSE],
                                    fml = fml_w, fam = fam, smoothing = smoothing, grid = grid,
                                    rng_kind = RNGkind(), opts = .refit_options())
      # each refit's warnings and messages are raised here as the loop below
      # raises them: in order, from the random number state that loop is in
      # when it fits that resample (so a handler that stops the call there
      # leaves the same state), and those of the fit inside the loop's catch
      # of a failed fit, so one that options(warn = 2) turns into an error
      # drops the refit as it does there
      say <- function(cnd) if (inherits(cnd, "warning")) warning(cnd) else message(cnd)
      moved <- FALSE
      for (b in seq_along(refits)) {
        r <- refits[[b]]
        if (!is.null(seeds[[b]]) && (length(r$fit_said) || length(r$predict_said) || !is.null(r$error))) {
          assign(".Random.seed", seeds[[b]], envir = env)
          moved <- TRUE
        }
        dropped <- tryCatch({
          for (cnd in r$fit_said) say(cnd)
          r$failed
        }, error = function(e) TRUE)
        if (dropped) next
        for (cnd in r$predict_said) say(cnd)
        if (!is.null(r$error)) stop(r$error)
        boot_fits[, b] <- r$curve
      }
      if (moved && !is.null(seeds[[n_boot]])) assign(".Random.seed", seeds[[n_boot]], envir = env)
    } else {
      for (b in seq_len(n_boot)) {
        fit_b <- fit_gam(df[sample.int(n, n, replace = TRUE), , drop = FALSE])
        if (is.null(fit_b)) next
        boot_fits[, b] <- fit_b$family$linkinv(
          as.numeric(stats::predict(fit_b, newdata = grid, type = "link"))
        )
      }
    }
    n_ok <- sum(!is.na(boot_fits[1, ]))
    if (n_ok < 2) {
      warning("Bootstrap ribbon failed (", n_ok, " usable resamples); using parametric interval")
      ci_method <- "parametric (Wald)"
      bootstrap <- FALSE
    } else {
      if (n_ok < n_boot) {
        warning(n_boot - n_ok, " of ", n_boot, " bootstrap resamples failed and were dropped")
      }
      grid$lwr <- apply(boot_fits, 1, stats::quantile, probs = 0.025, na.rm = TRUE)
      grid$upr <- apply(boot_fits, 1, stats::quantile, probs = 0.975, na.rm = TRUE)
    }
  }
  if (!bootstrap) {
    pr <- stats::predict(fit, newdata = grid, se.fit = TRUE, type = "link")
    grid$lwr <- linkinv(pr$fit - 1.96 * pr$se.fit)
    grid$upr <- linkinv(pr$fit + 1.96 * pr$se.fit)
  }

  # Ensure confidence bounds stay within [0,1] for binary fitness
  if (fitness_type == "binary") {
    grid$lwr <- pmax(grid$lwr, 0)
    grid$upr <- pmin(grid$upr, 1)
  }

  result <- list(
    model = fit,
    grid = grid,
    data = df[, c(trait_col, ".y")],
    trait = trait_col,
    fitness_type = fitness_type,
    family = family_name,
    dispersion = dispersion,
    k = k,
    basis = bs,
    smoothing = smoothing,
    spline_type = paste0(switch(bs, cr = "cubic regression spline", tp = "thin-plate spline", ps = "P-spline"),
                         " (", fit$method, ")"),
    ci_method = ci_method,
    n_obs = n_obs,
    fit_note = fit_note,
    group_used = group,
    trait_mean = z_mean,
    trait_sd = z_sd,
    surface_type = "correlated_fitness_univariate",
    note = "Univariate correlated fitness function (individual fitness)"
  )

  class(result) <- "univariate_fitness"

  return(result)
}
