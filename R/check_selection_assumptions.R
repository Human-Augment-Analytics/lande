# ======================================================
# check_selection_assumptions.R
# The checks a selection analysis rests on: multivariate normality of the
# traits (the condition for reading the gradients as the slope and curvature
# of the fitness surface),
# collinearity, rows per term, and the residuals or dispersion of the
# gradient models. Nothing is enforced.
# ======================================================

#' @noRd
# Mardia's (1970) multivariate skewness and kurtosis, with the small-sample
# factor for the skewness statistic. Rows are subsampled above `max_n`
# because the statistics need an n by n matrix; the subsample is drawn from
# a fixed seed, so the table is the same every run, and the caller's random
# numbers are left as they were.
.mardia <- function(X, max_n = 2000) {
  X <- as.matrix(X)
  X <- X[stats::complete.cases(X), , drop = FALSE]
  if (nrow(X) > max_n) X <- X[.with_seed(1, sample.int(nrow(X), max_n)), , drop = FALSE]
  n <- nrow(X)
  p <- ncol(X)
  Xc <- scale(X, scale = FALSE)
  S <- crossprod(Xc) / n
  Sinv <- tryCatch(solve(S), error = function(e) NULL)
  if (is.null(Sinv)) return(NULL)
  D <- Xc %*% Sinv %*% t(Xc)
  b1p <- sum(D^3) / n^2
  b2p <- sum(diag(D)^2) / n
  k <- (p + 1) * (n + 1) * (n + 3) / (n * ((n + 1) * (p + 1) - 6))
  skew <- n * k * b1p / 6
  p_skew <- stats::pchisq(skew, df = p * (p + 1) * (p + 2) / 6, lower.tail = FALSE)
  kurt <- (b2p - p * (p + 2)) / sqrt(8 * p * (p + 2) / n)
  p_kurt <- 2 * stats::pnorm(-abs(kurt))
  list(n = n, skewness = b1p, p_skewness = p_skew, kurtosis = b2p, p_kurtosis = p_kurt)
}

#' @noRd
# evaluate `expr` from a given seed and put the global random number state back
.with_seed <- function(seed, expr) {
  env <- globalenv()
  old <- if (exists(".Random.seed", envir = env, inherits = FALSE)) get(".Random.seed", envir = env) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = env) else assign(".Random.seed", old, envir = env))
  set.seed(seed)
  expr
}

#' @noRd
# Breusch-Pagan test of the squared residuals on the fitted values
.breusch_pagan <- function(fit) {
  r2 <- stats::residuals(fit)^2
  aux <- stats::lm(r2 ~ stats::fitted(fit))
  stat <- length(r2) * summary(aux)$r.squared
  c(statistic = stat, p = stats::pchisq(stat, df = 1, lower.tail = FALSE))
}

#' Check the assumptions behind a selection analysis
#'
#' Runs the checks a Lande and Arnold analysis rests on and returns them in
#' one table: multivariate normality of the traits, normality of each trait,
#' collinearity, rows per model term, and for the gradient models a residual
#' normality and heteroscedasticity test (continuous fitness), a separation
#' check (binary fitness) or the dispersion ratio (count fitness). Nothing is
#' enforced; the table is there to report alongside the gradients, as the
#' protocol of Palacio et al. (2019) asks.
#'
#' @details The regression gradients equal the covariance-adjusted
#'   differentials, \eqn{P^{-1} S}, by least squares whatever the trait
#'   distribution. Normal traits let them be read as the average slope and
#'   curvature of the fitness surface (Lande and Arnold 1983; Morrissey and
#'   Sakrejda 2013) and tie gamma to the change in the phenotypic covariance.
#'   Fitness doesn't need to be normal, and normality doesn't change how the
#'   gradients are computed.
#'   Mardia's (1970) skewness and kurtosis test multivariate normality, and
#'   each trait is also tested on its own with Shapiro and Wilk's test when
#'   there are 5000 rows or fewer.
#'
#'   In the package's simulation (validation/simulation.R in the GitHub
#'   repository), on a curved
#'   fitness surface, nominal 95\% intervals for beta covered 88\% of the time
#'   with normal traits, 67\% with log-normal ones and 77\% with symmetric
#'   heavy-tailed ones (92, 80 and 91\% by bootstrap). Most of that loss comes
#'   from the curvature, and HC3 standard errors
#'   (\code{selection_coefficients(se_type = "hc3")}) recover most of it.
#'   Estimating the SD of a heavy-tailed trait lowers coverage too, for gamma
#'   and, with skewed traits, for beta even on a straight surface (90\%), and
#'   HC3 doesn't help there. Skew also biases beta a little. Refitting with a
#'   transformed trait changes the scale selection is measured on, so report
#'   it with the untransformed result.
#'   Collinearity is the largest variance inflation factor from the linear
#'   gradient model. The rows-per-term check counts the individuals, or for
#'   binary fitness the rarer outcome, per term of the quadratic model, with
#'   ten as the working minimum. With more than 2000 individuals Mardia's
#'   statistics use a random subsample of 2000, drawn the same way every run.
#'   Like the gradient models, the logistic and count models fit an intercept
#'   per group, with unlabelled rows as one more group. A group in which every
#'   individual survived, or none did, is fitted at 0 or 1 by its intercept and
#'   is left out of the separation check. For counts the note names the model
#'   behind each set of p-values; the linear and quadratic models are checked
#'   for overdispersion separately. If the \code{performance} package is
#'   installed its heteroscedasticity and overdispersion tests are added next
#'   to the package's own, along with an R squared for the gradient model.
#'
#' @inheritParams selection_coefficients
#' @return A data frame of class \code{"selection_assumptions"} with one row
#'   per check: \code{check}, \code{statistic}, \code{p_value} and
#'   \code{note}.
#' @references Lande, R. and Arnold, S. J. (1983) The measurement of selection
#'   on correlated characters. Evolution 37, 1210-1226. Mardia, K. V. (1970)
#'   Measures of multivariate skewness and kurtosis with applications.
#'   Biometrika 57, 519-530. Morrissey, M. B. and Sakrejda, K. (2013)
#'   Unification of regression-based methods for the analysis of natural
#'   selection. Evolution 67, 2094-2100. Palacio, F. X., Ordano, M. and Benitez-Vieyra, S.
#'   (2019) Measuring natural selection on multivariate phenotypic traits: a
#'   protocol for verifiable and reproducible analyses of natural selection.
#'   Israel Journal of Ecology and Evolution 65, 130-136.
#' @export
#' @examples
#' check_selection_assumptions(bumpus, "survival", c("total_length", "weight", "humerus"))
check_selection_assumptions <- function(data,
                                        fitness_col,
                                        trait_cols,
                                        fitness_type = c("auto", "binary", "count", "continuous"),
                                        standardize = TRUE,
                                        group = NULL) {
  fitness_type <- match.arg(fitness_type)
  need <- c(fitness_col, trait_cols, group)
  absent <- setdiff(need, names(data))
  if (length(absent)) stop("Missing columns: ", paste(absent, collapse = ", "))

  # the rows the gradient models use: complete fitness and traits, with
  # unlabelled rows as one more group
  keep <- stats::complete.cases(data[, c(fitness_col, trait_cols), drop = FALSE])
  prep <- suppressWarnings(suppressMessages(prepare_selection_data(
    data[keep, , drop = FALSE], fitness_col, trait_cols, standardize = standardize, group = group,
    add_relative = TRUE, na_action = "none", name_relative = ".w"
  )))
  if (fitness_type == "auto") fitness_type <- detect_family(prep[[fitness_col]])$type
  if (!fitness_type %in% c("binary", "count")) fitness_type <- "continuous"
  # a group with zero mean fitness gets NA relative fitness; those rows are
  # outside the gradient models, so they are outside the counts here too
  prep <- prep[is.finite(prep$.w), , drop = FALSE]
  n <- nrow(prep)
  p <- length(trait_cols)
  rows <- list()
  add <- function(check, statistic = NA_real_, p_value = NA_real_, note = "") {
    rows[[length(rows) + 1]] <<- data.frame(check = check, statistic = statistic, p_value = p_value, note = note,
                                            stringsAsFactors = FALSE)
  }

  # traits
  X <- prep[, trait_cols, drop = FALSE]
  if (p >= 2) {
    m <- .mardia(X)
    if (!is.null(m)) {
      sub <- if (m$n < n) sprintf(" (random %d of %d rows)", m$n, n) else ""
      add("Multivariate normality: Mardia skewness", m$skewness, m$p_skewness,
          paste0(if (m$p_skewness < 0.05) "skewed" else "no evidence against normality", sub))
      add("Multivariate normality: Mardia kurtosis", m$kurtosis, m$p_kurtosis,
          if (m$p_kurtosis < 0.05) "tails differ from normal" else "no evidence against normality")
    }
  }
  for (t in trait_cols) {
    x <- X[[t]]
    if (n >= 3 && n <= 5000) {
      sw <- stats::shapiro.test(x)
      add(paste0("Normality of ", t, " (Shapiro-Wilk)"), unname(sw$statistic), sw$p.value,
          if (sw$p.value < 0.05) "not normal" else "")
    } else {
      z <- (x - mean(x)) / stats::sd(x)
      add(paste0("Skewness of ", t), mean(z^3), NA_real_, "Shapiro-Wilk needs 3 to 5000 rows")
    }
  }

  # collinearity
  if (p >= 2) {
    fml <- stats::as.formula(paste(".w ~", paste(trait_cols, collapse = " + ")))
    vif <- .compute_vif(stats::lm(fml, data = prep))
    if (!is.null(vif) && all(is.finite(vif))) {
      add("Collinearity: largest VIF", max(vif), NA_real_,
          if (max(vif) > 5) paste0("above 5 for ", trait_cols[which.max(vif)]) else "below 5")
    }
  }

  # rows per term of the quadratic model
  terms <- 2 * p + p * (p - 1) / 2
  events <- if (fitness_type == "binary") min(sum(prep[[fitness_col]] == 1), sum(prep[[fitness_col]] == 0)) else n
  per_term <- events / terms
  add(if (fitness_type == "binary") "Rarer outcome per quadratic term" else "Rows per quadratic term", per_term, NA_real_,
      if (per_term < 10) sprintf("%d %s for %d terms; treat gamma with care", events,
                                 if (fitness_type == "binary") "of the rarer outcome" else "rows", terms) else "")

  # the gradient models; the performance package adds its own tests and an R2 when installed
  has_perf <- requireNamespace("performance", quietly = TRUE)
  lin <- paste(trait_cols, collapse = " + ")
  # the logistic and count models take an intercept per group, as the gradient models do
  glm_data <- .add_glm_group(prep, group)
  glm_rhs <- if (".group" %in% names(glm_data)) paste(".group +", lin) else lin
  if (fitness_type == "continuous") {
    fit <- stats::lm(stats::as.formula(paste(".w ~", lin)), data = prep)
    r <- stats::residuals(fit)
    if (length(r) >= 3 && length(r) <= 5000) {
      sw <- stats::shapiro.test(r)
      add("Linear model residuals: normality (Shapiro-Wilk)", unname(sw$statistic), sw$p.value,
          if (sw$p.value < 0.05) "residuals not normal; bootstrap the intervals" else "")
    }
    bp <- .breusch_pagan(fit)
    add("Linear model residuals: heteroscedasticity (Breusch-Pagan)", bp[["statistic"]], bp[["p"]],
        if (bp[["p"]] < 0.05) "variance changes with the fitted value; bootstrap the intervals" else "")
    if (has_perf) {
      ph <- tryCatch(as.numeric(performance::check_heteroscedasticity(fit)), error = function(e) NA_real_)
      add("Linear model residuals: heteroscedasticity (performance)", NA_real_, ph,
          if (isTRUE(ph < 0.05)) "performance's test also finds the variance is not constant" else "")
    }
  } else if (fitness_type == "binary") {
    fit <- suppressWarnings(stats::glm(stats::as.formula(paste(fitness_col, "~", glm_rhs)), data = glm_data, family = stats::binomial()))
    coefs <- stats::coef(fit)[trait_cols]
    extreme <- max(abs(coefs), na.rm = TRUE)
    # a group in which everyone survived, or nobody did, is fitted at 0 or 1
    # by its intercept; only groups whose outcome varies can show separation
    y <- glm_data[[fitness_col]]
    varies <- if (".group" %in% names(glm_data)) {
      stats::ave(y, glm_data$.group, FUN = function(v) as.numeric(length(unique(v)) > 1)) == 1
    } else {
      rep(TRUE, length(y))
    }
    at_edge <- stats::fitted(fit) < 1e-6 | stats::fitted(fit) > 1 - 1e-6
    fitted_edge <- if (any(varies)) mean(at_edge[varies]) else 0
    add("Logistic model: separation", extreme, NA_real_,
        if (extreme > 10 || fitted_edge > 0.05) "coefficients or fitted values at the edge; possible complete separation" else "no sign of separation")
  } else {
    fit <- suppressWarnings(stats::glm(stats::as.formula(paste(fitness_col, "~", glm_rhs)), data = glm_data, family = stats::poisson()))
    disp <- sum(stats::residuals(fit, type = "pearson")^2) / stats::df.residual(fit)
    p_disp <- NA_real_
    if (has_perf) {
      od <- tryCatch(performance::check_overdispersion(fit), error = function(e) NULL)
      if (!is.null(od)) p_disp <- as.numeric(od$p_value)
    }
    # the linear and quadratic gradient models each move to a negative binomial
    # on their own dispersion, so name the model each set of gradients uses
    quad <- paste(c(trait_cols, paste0("I(", trait_cols, "^2)"),
                    if (p >= 2) utils::combn(trait_cols, 2, paste, collapse = ":")), collapse = " + ")
    family_of <- function(rhs) {
      rhs <- if (".group" %in% names(glm_data)) paste(".group +", rhs) else rhs
      tryCatch(attr(suppressWarnings(.fit_pvalue_glm(stats::as.formula(paste(fitness_col, "~", rhs)), glm_data, "count")),
                    "family_label"), error = function(e) NA_character_)
    }
    nb <- c(family_of(lin), family_of(quad)) %in% "negative binomial"
    note <- if (all(nb)) {
      "overdispersed; the p-values use a negative binomial model"
    } else if (nb[1]) {
      "overdispersed; the linear p-values use a negative binomial model, the quadratic ones a Poisson model"
    } else if (nb[2]) {
      "the quadratic model is overdispersed; its p-values use a negative binomial model, the linear ones a Poisson model"
    } else if (disp > 1.5) {
      "overdispersed, but the negative binomial fit failed; the p-values are from a Poisson model"
    } else {
      ""
    }
    add("Poisson model: dispersion ratio", disp, p_disp, note)
  }
  if (has_perf) {
    r2 <- tryCatch(performance::r2(fit), error = function(e) NULL)
    if (!is.null(r2) && length(r2)) {
      add(paste0("Model fit: ", names(r2)[1], " (performance)"), as.numeric(r2[[1]]), NA_real_, "")
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "fitness_type") <- fitness_type
  attr(out, "n") <- n
  class(out) <- c("selection_assumptions", "data.frame")
  out
}

#' Print the assumption checks
#'
#' @param x An object of class \code{"selection_assumptions"}.
#' @param ... Additional arguments (ignored).
#' @return The input object \code{x}, invisibly.
#' @export
print.selection_assumptions <- function(x, ...) {
  cat("Assumption checks (", attr(x, "fitness_type"), " fitness, n = ", attr(x, "n"), ")\n\n", sep = "")
  out <- data.frame(
    check = x$check,
    statistic = ifelse(is.na(x$statistic), "", formatC(x$statistic, digits = 3, format = "g")),
    p = ifelse(is.na(x$p_value), "", ifelse(x$p_value < 0.001, "< 0.001", formatC(x$p_value, digits = 3, format = "f"))),
    note = x$note,
    stringsAsFactors = FALSE
  )
  print(out, row.names = FALSE, right = FALSE)
  invisible(x)
}
