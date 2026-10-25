# ============================================================================
# extract_results.R
# Extract selection coefficients from model results
# ============================================================================

#' @noRd
# internal utility: a small table as a message, so it can be silenced with
# suppressMessages() like every other note the functions print
.msg_table <- function(title, x) {
  message(title, "\n", paste(utils::capture.output(print(x)), collapse = "\n"))
}

#' @noRd
# internal utility: variance inflation factors, with the collinearity warning
.compute_vif <- function(fit) {
  if (!requireNamespace("car", quietly = TRUE)) {
    return(NULL)
  }
  # VIF is only defined with two or more predictors; a single-trait model has
  # nothing to inflate, so return NULL rather than warn.
  if (length(attr(stats::terms(fit), "term.labels")) < 2) {
    return(NULL)
  }
  vif_vals <- tryCatch(
    car::vif(fit),
    error = function(e) {
      warning("VIF calculation failed: ", e$message)
      NULL
    }
  )
  # an exact fit (no residual degrees of freedom) has NaN VIFs, which say
  # nothing about collinearity
  if (!is.null(vif_vals) && any(vif_vals > 5, na.rm = TRUE)) {
    warning("Collinear traits (VIF above 5) may inflate the standard errors")
  }
  vif_vals
}

#' @noRd
# internal utility: warn when the design is under-powered for its parameters
.warn_small_sample <- function(n, n_params) {
  if (n < n_params * 2) {
    warning(
      "Sample size (", n, ") may be too small for ", n_params,
      " parameters"
    )
  }
}

#' @noRd
# internal utility: reference level used when predicting from a model that
# includes a group term. Numeric groups use the level closest to the median of
# the observed levels; factor/character groups use the most common level.
# Missing labels are ignored so a single NA cannot leave the choice empty.
.reference_group <- function(x) {
  x <- x[!is.na(as.character(x))]
  if (!length(x)) {
    stop("Group column has no non-missing values")
  }
  if (is.numeric(x)) {
    lv <- unique(x)
    lv[which.min(abs(lv - stats::median(lv)))]
  } else {
    names(sort(table(x), decreasing = TRUE))[1]
  }
}

#' @noRd
# internal utility: does a raw fitness column look like the type claimed?
.is_raw_fitness <- function(raw, fitness_type) {
  if (fitness_type == "binary") {
    return(all(raw %in% c(0, 1)))
  }
  all(raw >= 0 & abs(raw - round(raw)) < 1e-8)
}

#' @noRd
# internal utility: the mgcv family for count fitness in the spline and the
# surface.
.count_family <- function(count_family) {
  switch(count_family,
         poisson = stats::poisson("log"),
         quasipoisson = stats::quasipoisson("log"),
         nb = mgcv::nb())
}

#' @noRd
# internal utility: the Pearson dispersion of a count GAM. A Poisson fit warns
# above 1.5, the rule the gradient p-values use.
.check_dispersion <- function(fit, count_family) {
  disp <- sum(stats::residuals(fit, type = "pearson")^2) / fit$df.residual
  if (count_family == "poisson" && is.finite(disp) && disp > 1.5) {
    warning("Counts are overdispersed (dispersion ", round(disp, 2),
            "); the Poisson standard errors are too small, and count_family = \"quasipoisson\" corrects them")
  }
  disp
}

#' @noRd
# internal utility: the group as a factor for the p-value GLM when there are
# two or more; unlabelled rows count as one more group, as in the prepared
# data.
.add_glm_group <- function(data, group) {
  if (is.null(group) || !group %in% names(data)) return(data)
  g <- droplevels(addNA(factor(data[[group]]), ifany = TRUE))
  if (nlevels(g) > 1) data$.group <- g
  data
}

#' @noRd
# internal utility: the standard errors a selection model reports, "ols" or
# "hc3", in either case
.se_type_arg <- function(se_type) {
  match.arg(tolower(se_type[1]), c("ols", "hc3"))
}

#' @noRd
# internal utility: swap the least-squares standard errors in a summary.lm for
# heteroscedasticity-consistent HC3 ones (MacKinnon and White 1985) and redo
# the t tests on them, with the residual degrees of freedom. HC3 divides each
# residual by one minus its leverage, so it is undefined at leverage 1.
.hc3_summary <- function(fit, sm) {
  keep <- !is.na(stats::coef(fit))
  X <- stats::model.matrix(fit)[, keep, drop = FALSE]
  h <- stats::hatvalues(fit)
  if (any(h > 1 - 1e-8)) {
    warning("HC3 standard errors are undefined when a row has leverage 1")
  }
  e <- stats::residuals(fit) / (1 - h)
  B <- solve(crossprod(X))
  se <- sqrt(diag(B %*% crossprod(X * e) %*% B))
  cf <- sm$coefficients
  cf[, "Std. Error"] <- se[rownames(cf)]
  cf[, "t value"] <- cf[, "Estimate"] / cf[, "Std. Error"]
  cf[, "Pr(>|t|)"] <- 2 * stats::pt(-abs(cf[, "t value"]), df = fit$df.residual)
  sm$coefficients <- cf
  sm
}

#' @noRd
# internal utility: the GLM that supplies p-values when OLS residuals cannot
# be trusted. Binary fitness gets a logistic model. Counts get a Poisson
# model, swapped for a negative binomial when the Pearson dispersion is above
# 1.5, since overdispersed counts make Poisson p-values too small.
.fit_pvalue_glm <- function(formula, data, fitness_type) {
  if (fitness_type == "binary") {
    fit <- stats::glm(formula, data = data, family = stats::binomial())
    attr(fit, "family_label") <- "binomial(logit)"
    return(fit)
  }
  fit <- stats::glm(formula, data = data, family = stats::poisson())
  disp <- sum(stats::residuals(fit, type = "pearson")^2) / fit$df.residual
  label <- "poisson(log)"
  if (is.finite(disp) && disp > 1.5) {
    nb <- tryCatch(
      suppressWarnings(MASS::glm.nb(formula, data = data)),
      error = function(e) NULL
    )
    if (is.null(nb)) {
      warning(
        "Counts are overdispersed (dispersion ", round(disp, 2),
        ") but the negative binomial fit failed; Poisson p-values kept"
      )
    } else {
      fit <- nb
      label <- "negative binomial"
    }
  }
  attr(fit, "dispersion") <- disp
  attr(fit, "family_label") <- label
  fit
}

#' @noRd
# internal utility: get coefficient-level p-value column name from summary()
.p_col_from_summary <- function(coef_mat) {
  pcols <- intersect(colnames(coef_mat), c("Pr(>|t|)", "Pr(>|z|)"))
  if (length(pcols)) pcols[1] else NA_character_
}

#' @noRd
# internal utility: safely extract ANOVA p-value (car::Anova)
.ano_p <- function(anova_obj, term) {
  if (is.null(anova_obj)) {
    return(NA_real_)
  }
  rn <- rownames(anova_obj)
  if (is.null(rn) || !term %in% rn) {
    return(NA_real_)
  }
  pcols <- intersect(colnames(anova_obj), c("Pr(>F)", "Pr(>Chisq)"))
  if (!length(pcols)) {
    return(NA_real_)
  }
  as.numeric(anova_obj[term, pcols[1]])
}
#' @noRd
# Helper to get appropriate summary and p-value column
.get_summary_and_pcol <- function(results) {
  # Binary and count fitness carry a second model whose p-values replace the
  # OLS ones; `is_binary` below means "p-values come from that GLM".
  if (results$fitness_type %in% c("binary", "count")) {
    # Check if summary$glm exists
    if (is.null(results$summary$glm)) {
      sm <- results$summary$ols
      pcol <- .p_col_from_summary(coef(sm))
      if (is.na(pcol)) pcol <- "Pr(>|t|)"

      return(list(
        coef_mat = coef(sm),
        pcol = pcol,
        anova = results$anova,
        is_binary = FALSE
      ))
    }

    sm_ols <- results$summary$ols
    sm_glm <- results$summary$glm

    coef_ols <- coef(sm_ols)
    coef_glm <- coef(sm_glm)

    if (is.null(coef_glm)) {
      stop("GLM coefficients are NULL. GLM model may have failed.")
    }

    pcol_glm <- .p_col_from_summary(coef(sm_glm))
    if (is.na(pcol_glm)) pcol_glm <- "Pr(>|z|)"

    if (!is.matrix(coef_glm)) {
      if (is.vector(coef_glm) && length(coef_glm) > 0) {
        coef_glm <- matrix(coef_glm, nrow = 1)
        rownames(coef_glm) <- names(coef_glm)
      }
    }

    return(list(
      coef_mat_ols = coef_ols,
      coef_mat_glm = coef_glm,
      pcol_glm = pcol_glm,
      anova = results$anova,
      is_binary = TRUE
    ))
  } else {
    # Continuous case
    sm <- results$summary
    pcol <- .p_col_from_summary(coef(sm))
    if (is.na(pcol)) pcol <- "Pr(>|t|)"

    return(list(
      coef_mat = coef(sm),
      pcol = pcol,
      anova = results$anova,
      is_binary = FALSE
    ))
  }
}

#' Extract linear selection coefficients
#'
#' @param trait_cols A character vector of trait column names.
#' @param results A model results object returned by \code{analyze_linear_selection()}.
#'
#' @return A data frame with linear selection coefficients and statistics.
#' @examples
#' prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
#' fit <- analyze_linear_selection(prep, "survival", c("total_length", "weight"), "binary")
#' extract_linear_coefficients(c("total_length", "weight"), fit)
#' @export
extract_linear_coefficients <- function(trait_cols, results) {
  obj <- .get_summary_and_pcol(results)

  if (obj$is_binary) {
    # Binary case: gradients from OLS, p-values from GLM
    coef_ols <- obj$coef_mat_ols
    coef_glm <- obj$coef_mat_glm
    pcol_glm <- obj$pcol_glm
    anova_obj <- obj$anova

    keep <- intersect(rownames(coef_ols), trait_cols)

    rows <- lapply(keep, function(t) {
      # p-value from GLM
      if (t %in% rownames(coef_glm)) {
        p_val <- coef_glm[t, pcol_glm]
      } else {
        p_val <- .ano_p(anova_obj, t)
      }

      data.frame(
        Term = t,
        Type = "Linear",
        Beta_Coefficient = coef_ols[t, "Estimate"],
        Standard_Error = coef_ols[t, "Std. Error"],
        P_Value = as.numeric(p_val),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })

    if (!length(rows)) {
      return(data.frame(
        Term = character(), Type = character(),
        Beta_Coefficient = numeric(), Standard_Error = numeric(),
        P_Value = numeric(), check.names = FALSE
      ))
    }
    return(do.call(rbind, rows))
  } else {
    # Continuous case: everything from OLS
    coef_mat <- obj$coef_mat
    pcol <- obj$pcol
    anova_obj <- obj$anova

    keep <- intersect(rownames(coef_mat), trait_cols)

    rows <- lapply(keep, function(t) {
      p_val <- if (!is.na(pcol) && t %in% rownames(coef_mat)) {
        coef_mat[t, pcol]
      } else {
        .ano_p(anova_obj, t)
      }

      data.frame(
        Term = t,
        Type = "Linear",
        Beta_Coefficient = coef_mat[t, "Estimate"],
        Standard_Error = coef_mat[t, "Std. Error"],
        P_Value = as.numeric(p_val),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })

    if (!length(rows)) {
      return(data.frame(
        Term = character(), Type = character(),
        Beta_Coefficient = numeric(), Standard_Error = numeric(),
        P_Value = numeric(), check.names = FALSE
      ))
    }
    return(do.call(rbind, rows))
  }
}

#' Extract quadratic selection coefficients
#'
#' @param trait_cols A character vector of trait column names.
#' @param results A model results object returned by \code{analyze_nonlinear_selection()}.
#'
#' @return A data frame with quadratic selection coefficients (doubled estimates) and statistics.
#' @examples
#' prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
#' fit <- analyze_nonlinear_selection(prep, "survival", c("total_length", "weight"), "binary")
#' extract_quadratic_coefficients(c("total_length", "weight"), fit)
#' @export
extract_quadratic_coefficients <- function(trait_cols, results) {
  obj <- .get_summary_and_pcol(results)

  if (obj$is_binary) {
    # Binary case: gradients from OLS, p-values from GLM
    coef_ols <- obj$coef_mat_ols
    coef_glm <- obj$coef_mat_glm
    pcol_glm <- obj$pcol_glm
    anova_obj <- obj$anova

    rows <- lapply(trait_cols, function(t) {
      term <- paste0("I(", t, "^2)")
      if (!term %in% rownames(coef_ols)) {
        return(NULL)
      }

      # Multiply estimate and SE by 2 for gamma_ii
      est <- 2 * coef_ols[term, "Estimate"]
      se <- 2 * coef_ols[term, "Std. Error"]

      # p-value from GLM
      if (term %in% rownames(coef_glm)) {
        p_val <- coef_glm[term, pcol_glm]
      } else {
        p_val <- .ano_p(anova_obj, term)
      }

      data.frame(
        Term = paste0(t, "\u00b2"),
        Type = "Quadratic",
        Beta_Coefficient = as.numeric(est),
        Standard_Error = as.numeric(se),
        P_Value = as.numeric(p_val),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })

    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) {
      return(data.frame(
        Term = character(), Type = character(),
        Beta_Coefficient = numeric(), Standard_Error = numeric(),
        P_Value = numeric(), check.names = FALSE
      ))
    }
    return(do.call(rbind, rows))
  } else {
    # Continuous case: everything from OLS
    coef_mat <- obj$coef_mat
    pcol <- obj$pcol
    anova_obj <- obj$anova

    rows <- lapply(trait_cols, function(t) {
      term <- paste0("I(", t, "^2)")
      if (!term %in% rownames(coef_mat)) {
        return(NULL)
      }

      est <- 2 * coef_mat[term, "Estimate"]
      se <- 2 * coef_mat[term, "Std. Error"]
      p_val <- if (!is.na(pcol) && term %in% rownames(coef_mat)) {
        coef_mat[term, pcol]
      } else {
        .ano_p(anova_obj, term)
      }

      data.frame(
        Term = paste0(t, "\u00b2"),
        Type = "Quadratic",
        Beta_Coefficient = as.numeric(est),
        Standard_Error = as.numeric(se),
        P_Value = as.numeric(p_val),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })

    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) {
      return(data.frame(
        Term = character(), Type = character(),
        Beta_Coefficient = numeric(), Standard_Error = numeric(),
        P_Value = numeric(), check.names = FALSE
      ))
    }
    return(do.call(rbind, rows))
  }
}


#' Extract correlational selection coefficients
#'
#' @param trait_cols A character vector of trait column names.
#' @param results A model results object returned by \code{analyze_nonlinear_selection()}.
#'
#' @return A data frame with correlational (interaction) selection coefficients and statistics.
#' @examples
#' prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
#' fit <- analyze_nonlinear_selection(prep, "survival", c("total_length", "weight"), "binary")
#' extract_interaction_coefficients(c("total_length", "weight"), fit)
#' @export
extract_interaction_coefficients <- function(trait_cols, results) {
  if (length(trait_cols) < 2) {
    return(data.frame(
      Term = character(), Type = character(),
      Beta_Coefficient = numeric(), Standard_Error = numeric(),
      P_Value = numeric(), check.names = FALSE
    ))
  }

  obj <- .get_summary_and_pcol(results)
  pairs <- utils::combn(trait_cols, 2, simplify = FALSE)

  if (obj$is_binary) {
    # Binary case: gradients from OLS, p-values from GLM
    coef_ols <- obj$coef_mat_ols
    coef_glm <- obj$coef_mat_glm
    pcol_glm <- obj$pcol_glm
    anova_obj <- obj$anova

    rows <- lapply(pairs, function(p) {
      raw_term <- paste0(p[1], ":", p[2])
      if (!raw_term %in% rownames(coef_ols)) {
        return(NULL)
      }

      est <- coef_ols[raw_term, "Estimate"]
      se <- coef_ols[raw_term, "Std. Error"]

      # p-value from GLM
      if (raw_term %in% rownames(coef_glm)) {
        p_val <- coef_glm[raw_term, pcol_glm]
      } else {
        p_val <- .ano_p(anova_obj, raw_term)
      }

      data.frame(
        Term = paste(p[1], "\u00d7", p[2]),
        Type = "Correlational",
        Beta_Coefficient = as.numeric(est),
        Standard_Error = as.numeric(se),
        P_Value = as.numeric(p_val),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })

    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) {
      return(data.frame(
        Term = character(), Type = character(),
        Beta_Coefficient = numeric(), Standard_Error = numeric(),
        P_Value = numeric(), check.names = FALSE
      ))
    }
    return(do.call(rbind, rows))
  } else {
    # Continuous case: everything from OLS
    coef_mat <- obj$coef_mat
    pcol <- obj$pcol
    anova_obj <- obj$anova

    rows <- lapply(pairs, function(p) {
      raw_term <- paste0(p[1], ":", p[2])
      if (!raw_term %in% rownames(coef_mat)) {
        return(NULL)
      }

      est <- coef_mat[raw_term, "Estimate"]
      se <- coef_mat[raw_term, "Std. Error"]
      p_val <- if (!is.na(pcol) && raw_term %in% rownames(coef_mat)) {
        coef_mat[raw_term, pcol]
      } else {
        .ano_p(anova_obj, raw_term)
      }

      data.frame(
        Term = paste(p[1], "\u00d7", p[2]),
        Type = "Correlational",
        Beta_Coefficient = as.numeric(est),
        Standard_Error = as.numeric(se),
        P_Value = as.numeric(p_val),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })

    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) {
      return(data.frame(
        Term = character(), Type = character(),
        Beta_Coefficient = numeric(), Standard_Error = numeric(),
        P_Value = numeric(), check.names = FALSE
      ))
    }
    return(do.call(rbind, rows))
  }
}
