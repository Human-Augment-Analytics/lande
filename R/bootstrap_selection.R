# ============================================================================
# bootstrap_selection
#
# Bootstrap standard errors and confidence intervals for selection gradients.
# Individuals are resampled with replacement and the whole estimation procedure
# (standardisation, relativisation, OLS) is repeated on each resample, so the
# intervals include standardising and relativising as well as the fit. Only
# the least-squares gradients are refitted; the p-value models are fitted
# once, on the full data.
# ============================================================================

#' Bootstrap selection gradients
#'
#' Resamples individuals with replacement and re-estimates the selection
#' gradients on each resample to obtain bootstrap standard errors and percentile
#' confidence intervals for every coefficient (beta, gamma, gamma_ij).
#'
#' @inheritParams selection_coefficients
#' @param n_boot Integer number of bootstrap resamples. Default is 1000.
#' @param conf Confidence level for the percentile interval. Default is 0.95.
#'
#' @details Each resample goes through the whole procedure again: with
#'   \code{standardize = TRUE} the traits are restandardised, and fitness is made
#'   relative to the resample's own mean, before both least-squares models are
#'   refitted, so the intervals carry the uncertainty in the standardisation as
#'   well as in the fit. Only those two models are refitted on each resample;
#'   the GLMs that supply the p-values (logistic for binary fitness, Poisson or
#'   negative binomial for counts) are fitted once, on the full data, so a GLM
#'   that would fail on some resample cannot remove it. Resamples whose
#'   least-squares fit fails are dropped, and \code{N_Boot} counts the ones
#'   that were used. If the column names keep this shortcut from reproducing
#'   the point estimate exactly (a trait named like an R expression, say),
#'   each resample is refitted with all of \code{selection_coefficients()}
#'   instead, which is slower and drops a resample whose GLM fails.
#'
#'   When \code{group} is given, individuals are resampled within each
#'   group so every resample keeps the original group sizes. A resample in
#'   which a trait has no variance within a group (so it could only be centred)
#'   is discarded rather than fitted on a degenerate value. The point estimate
#'   is the ordinary \code{selection_coefficients()} fit on the full data and
#'   its warnings are reported as usual.
#'
#'   With continuous fitness the intervals for beta can be too narrow when the
#'   traits have heavy tails (skewed or not) and the linear model leaves out
#'   curvature, though less so than the parametric ones (see
#'   \code{check_selection_assumptions()}). For survival and counts in the
#'   package's simulation they covered a little less than the parametric
#'   ones, 92 to 93\% for beta. Refitting
#'   on a transformed trait changes the scale of the gradients and is only a
#'   check.
#'
#' @return A data frame with one row per coefficient and columns \code{Term},
#'   \code{Type}, \code{Estimate} (point estimate on the full data),
#'   \code{Boot_SE}, \code{CI_lower}, \code{CI_upper}, \code{P_Value}, and
#'   \code{N_Boot} (usable resamples for that coefficient). \code{P_Value} is
#'   the parametric p-value from the point-estimate fit (the OLS t-test for
#'   continuous fitness, the Wald test of the logistic GLM for binary fitness
#'   or of the Poisson or negative binomial GLM for counts), not a bootstrap
#'   p-value.
#' @export
#'
#' @examples
#' set.seed(1)
#' bootstrap_selection(bumpus, "survival", c("total_length", "weight"), n_boot = 50)
bootstrap_selection <- function(data,
                                fitness_col,
                                trait_cols,
                                fitness_type = c("auto", "binary", "count", "continuous"),
                                standardize = TRUE,
                                group = NULL,
                                use_relative_for_fit = TRUE,
                                n_boot = 1000,
                                conf = 0.95) {
  fitness_type <- match.arg(fitness_type)
  if (!is.numeric(conf) || length(conf) != 1L || is.na(conf) || conf <= 0 || conf >= 1) {
    stop("`conf` must be a single number between 0 and 1")
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1L || is.na(n_boot) ||
      n_boot < 2 || n_boot != round(n_boot)) {
    stop("`n_boot` must be a single whole number of at least 2")
  }
  n_boot <- as.integer(n_boot)
  if (!is.null(group) && !group %in% names(data)) {
    stop("Group column '", group, "' not found in data")
  }

  fit_once <- function(d) {
    selection_coefficients(
      d, fitness_col, trait_cols,
      fitness_type = fitness_type,
      standardize = standardize,
      group = group,
      use_relative_for_fit = use_relative_for_fit
    )
  }

  # Point estimate on the full data, with its p-values; its warnings are shown.
  point <- suppressMessages(fit_once(data))
  keys <- paste(point$Term, point$Type)
  n <- nrow(data)

  # The shortcut: a resample refits only the two least-squares models, with
  # .ols_gradients() on the raw columns of the resampled rows, and none of the
  # GLMs and checks that give the point estimate its p-values. NULL for an
  # unusable resample: a trait lost its variance in some group (preparing the
  # data then centres only, which would put a false z = 0 into the fit), or
  # the fit failed.
  src <- if (is.data.frame(data)) data else data.frame(data, check.names = FALSE)
  cols <- lapply(stats::setNames(nm = unique(c(fitness_col, trait_cols))), function(v) src[[v]])
  codes <- .group_codes(src, group)
  pos <- match(keys, .gradient_keys(trait_cols))
  gradients <- function(rows, drop_degenerate) {
    g <- tryCatch(
      suppressWarnings(.ols_gradients(
        lapply(cols, `[`, rows), fitness_col, trait_cols, codes[rows],
        standardize = standardize,
        use_relative_for_fit = use_relative_for_fit,
        drop_degenerate = drop_degenerate
      )),
      error = function(e) NULL
    )
    if (is.null(g)) NULL else as.vector(g[pos])
  }

  # The shortcut is used only if it gives the point estimate's gradients bit
  # for bit on the full data. Column names can keep it from doing so: the
  # model formula reads a trait named "log(x)" as log() of x, and the
  # preparation overwrites a trait named like the relative-fitness column.
  # Then every resample goes through selection_coefficients(), as it always
  # did; this is also the check that the two never drift apart.
  fit_resample <- if (identical(gradients(seq_len(n), FALSE), as.vector(point$Beta_Coefficient))) {
    function(rows) gradients(rows, drop_degenerate = TRUE)
  } else {
    function(rows) {
      d <- data[rows, , drop = FALSE]
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
      if (degenerate || is.null(res)) return(NULL)
      res$Beta_Coefficient[match(keys, paste(res$Term, res$Type))]
    }
  }

  # Resample within groups when there are groups, so each resample keeps the
  # original group sizes. Rows with no group label are their own group, as
  # when the data are prepared; split() would drop them, even after addNA(),
  # so it splits on the level codes.
  strata <- if (is.null(group)) rep(1L, n) else codes
  strata_rows <- split(seq_len(n), strata, drop = TRUE)
  resample_rows <- function() {
    unlist(lapply(strata_rows, function(r) r[sample.int(length(r), length(r), replace = TRUE)]),
           use.names = FALSE)
  }

  boot <- matrix(NA_real_, nrow = length(keys), ncol = n_boot,
                 dimnames = list(keys, NULL))

  for (b in seq_len(n_boot)) {
    res <- fit_resample(resample_rows())
    if (is.null(res)) next
    boot[, b] <- res
  }

  # Usable resamples per coefficient. A term that drops out of some fits (e.g. a
  # collinear interaction) has fewer usable resamples than the rest, so count
  # each row on its own instead of reading only the first.
  n_used <- rowSums(!is.na(boot))
  n_min <- min(n_used)
  if (n_min < 2) {
    stop("Bootstrap failed: fewer than 2 usable resamples for at least one term")
  }
  if (n_min < n_boot) {
    warning(n_boot - n_min, " of ", n_boot,
            " resamples were unusable for at least one term and were dropped")
  }

  a <- (1 - conf) / 2
  out <- data.frame(
    Term = point$Term,
    Type = point$Type,
    Estimate = point$Beta_Coefficient,
    Boot_SE = apply(boot, 1, stats::sd, na.rm = TRUE),
    CI_lower = apply(boot, 1, stats::quantile, probs = a, na.rm = TRUE),
    CI_upper = apply(boot, 1, stats::quantile, probs = 1 - a, na.rm = TRUE),
    P_Value = point$P_Value,
    N_Boot = as.integer(n_used),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(out) <- NULL
  attr(out, "n_boot") <- n_min
  attr(out, "conf") <- conf
  out
}

#' @noRd
# internal utility: "Term Type" for every gradient of the second-order model,
# named as selection_coefficients() names them and in the order
# .ols_gradients() returns them: linear, quadratic, then each pair of traits
# in combn() order.
.gradient_keys <- function(trait_cols) {
  pairs <- if (length(trait_cols) >= 2) utils::combn(trait_cols, 2, simplify = FALSE) else list()
  c(paste(trait_cols, "Linear"),
    paste(paste0(trait_cols, "\u00b2"), "Quadratic"),
    vapply(pairs, function(p) paste(paste(p[1], "\u00d7", p[2]), "Correlational"), character(1)))
}

#' @noRd
# internal utility: each row's group as an integer code, with a missing label
# as one more group, as prepare_selection_data() groups the rows; NULL
# without a group.
.group_codes <- function(data, group) {
  if (is.null(group)) return(NULL)
  as.integer(addNA(factor(data[[group]]), ifany = TRUE))
}

#' @noRd
# internal utility: the selection gradients alone, by least squares, with
# none of the p-value models; bootstrap_selection() refits them on every
# resample.
#
# `cols` is a named list of the raw fitness and trait columns and `codes` the
# group of each row (from .group_codes()). The traits are standardised and
# fitness made relative within groups, over the rows with complete fitness
# and traits, with prepare_selection_data()'s own utilities (.scaling_within()
# and .standardise(), .mean_fitness_within() and .relativise()); then the
# rows the gradient models keep are fitted, linear and quadratic, on the
# design matrices lm() builds from the formulas
# in analyze_linear_selection() and analyze_nonlinear_selection(), with the
# same QR routine and tolerance, and the quadratic coefficients are doubled.
#
# Returns the gradients in the order of .gradient_keys(trait_cols), NA for a
# term dropped as aliased, with attribute "degenerate" TRUE when a trait has
# no variance in some group (prepare_selection_data() warns and only centres
# it there). With drop_degenerate = TRUE such data return NULL instead. NULL
# also when fitting relative fitness and no group has a usable mean fitness.
.ols_gradients <- function(cols, fitness_col, trait_cols, codes = NULL,
                           standardize = TRUE, use_relative_for_fit = TRUE,
                           drop_degenerate = FALSE) {
  w <- cols[[fitness_col]]
  if (is.logical(w) || is.integer(w)) cols[[fitness_col]] <- as.numeric(w)
  cc <- !is.na(cols[[fitness_col]])
  for (t in trait_cols) cc <- cc & !is.na(cols[[t]])
  groups <- .group_rows(codes, length(cc))

  degenerate <- FALSE
  if (standardize) {
    for (t in trait_cols) {
      s <- .scaling_within(cols[[t]], cc, groups)
      unscalable <- !any(s$scalable[cc])
      if (unscalable || any(!s$scalable[cc])) {
        if (drop_degenerate) return(NULL)
        degenerate <- TRUE
        if (unscalable) next
      }
      cols[[t]] <- .standardise(cols[[t]], s)
    }
  }

  if (use_relative_for_fit) {
    m <- .mean_fitness_within(cols[[fitness_col]], cc, groups)
    if (!any(m$usable[cc])) return(NULL)
    y <- .relativise(cols[[fitness_col]], m)
  } else {
    y <- cols[[fitness_col]]
  }

  # the rows the gradient models keep: complete response, fitness and traits
  keep <- !is.na(y) & !is.na(cols[[fitness_col]])
  for (t in trait_cols) keep <- keep & !is.na(cols[[t]])
  y <- as.double(y[keep])
  if (!length(y)) stop("0 (non-NA) cases")
  # doubles before any product, so integer traits cannot overflow
  z <- lapply(cols[trait_cols], function(v) as.double(v[keep]))
  k <- length(z)
  one <- rep(1, length(y))
  # model.matrix() columns: intercept, traits, I(z^2) (R's z^2 is z * z), and
  # the products of each pair in the order of the formula's interactions
  sq <- lapply(z, function(v) v^2)
  cross <- if (k >= 2) {
    utils::combn(k, 2, function(p) z[[p[1]]] * z[[p[2]]], simplify = FALSE)
  } else {
    list()
  }
  lin <- .ls_coef(c(list(one), z), y)
  quad <- .ls_coef(c(list(one), z, sq, cross), y)
  out <- c(lin[1L + seq_len(k)], 2 * quad[1L + k + seq_len(k)], quad[-seq_len(1L + 2L * k)])
  attr(out, "degenerate") <- degenerate
  out
}

#' @noRd
# internal utility: least-squares coefficients of y on the given columns,
# with lm.fit()'s QR (same routine and tolerance) and its NA for aliased
# columns, without building an lm object.
.ls_coef <- function(columns, y) {
  X <- matrix(unlist(columns, use.names = FALSE), nrow = length(y))
  fit <- stats::.lm.fit(X, y)
  cf <- fit$coefficients
  p <- ncol(X)
  if (fit$rank < p) cf[(fit$rank + 1L):p] <- NA
  if (fit$pivoted) cf[fit$pivot] <- cf
  cf
}
