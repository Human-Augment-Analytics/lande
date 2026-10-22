# ============================================================================
# prepare_selection_data
#
# Purpose:
#   1. Check fitness and traits columns
#   2. Handle missing values (including group column)
#   3. Standardize trait variables (optionally within groups)
#   4. Compute relative fitness (optionally within groups)
#
# Groups:
#   For studies spanning several years, sites or populations, standardise the
#   traits within each group (year, site) rather than over all individuals
#   pooled, so each individual is compared with others in the same year and
#   place. Relative fitness is computed within the same group.
#
#   This approach is equivalent to including group fixed effects in the
#   regression models. However, it does NOT allow testing for group-by-trait
#   interactions (i.e., whether selection varies across years). For that,
#   researchers should fit separate models per group or include interaction
#   terms in the regression.
#
# Parameters:
#   data           : data frame containing fitness and trait measurements
#   fitness_col    : name of the fitness column (character)
#   trait_cols     : vector of trait column names
#   standardize    : logical; if TRUE, standardize traits to mean 0, SD 1
#   group          : optional grouping variable (e.g., "year", "site", "population")
#                    When specified, standardization and relative fitness are
#                    calculated separately within each group level.
#   add_relative   : logical; if TRUE, add a relative fitness column
#   na_action      : how to handle NAs ("warn", "drop", "none")
#   name_relative  : name for the relative fitness column
#
# Returns:
#   Modified data frame with:
#     - Standardized traits (if standardize = TRUE)
#     - Relative fitness column (if add_relative = TRUE)
#     - Potentially dropped NA rows (if na_action = "drop")
# ============================================================================

#' Prepare data for selection analysis
#'
#' Cleans, standardizes, and calculates relative fitness for trait and fitness data.
#' Standardizations and relative fitness calculations can be performed within groups.
#'
#' @param data A data frame containing fitness and trait measurements.
#' @param fitness_col A string specifying the name of the fitness column.
#' @param trait_cols A character vector of trait column names.
#' @param standardize Logical indicating whether to standardize traits to mean 0 and SD 1. Default is \code{TRUE}.
#' @param group Optional string specifying a grouping variable.
#' @param add_relative Logical indicating whether to add a relative fitness column. Default is \code{TRUE}.
#' @param na_action A string specifying how to handle missing values: \code{"warn"}, \code{"drop"}, or \code{"none"}.
#' @param name_relative A string specifying the name for the relative fitness column. Default is \code{"relative_fitness"}.
#'
#' @return A modified data frame ready for selection analysis.
#' @examples
#' prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
#' head(prep[, c("survival", "relative_fitness", "total_length", "weight")])
#'
#' # standardised and relativised within each sex
#' by_sex <- prepare_selection_data(bumpus, "survival", "total_length", group = "sex")
#' tapply(by_sex$total_length, by_sex$sex, mean)
#' @export
prepare_selection_data <- function(data,
                                   fitness_col,
                                   trait_cols,
                                   standardize = TRUE,
                                   group = NULL,
                                   add_relative = TRUE,
                                   na_action = c("warn", "drop", "none"),
                                   name_relative = "relative_fitness") {
  # Input validation
  na_action <- match.arg(na_action)

  df <- data.frame(data, check.names = FALSE)

  if (!is.character(fitness_col) || length(fitness_col) != 1L) {
    stop("`fitness_col` must be a single column name (character)")
  }

  if (!all(c(fitness_col, trait_cols) %in% names(df))) {
    stop("`fitness_col` or some `trait_cols` not found in data")
  }

  # Coerce fitness to numeric if it's logical or integer
  if (is.logical(df[[fitness_col]]) || is.integer(df[[fitness_col]])) {
    df[[fitness_col]] <- as.numeric(df[[fitness_col]])
  }

  # Check if group column exists
  if (!is.null(group)) {
    if (!group %in% names(df)) {
      stop("Group column '", group, "' not found in data")
    }
    message("Standardising and computing relative fitness within groups: '", group, "'")
  }

  # Handle missing values
  cols_check <- c(fitness_col, trait_cols)
  if (!is.null(group)) {
    cols_check <- c(cols_check, group)
  }

  na_rows <- !stats::complete.cases(df[, cols_check, drop = FALSE])

  if (any(na_rows)) {
    n_bad <- sum(na_rows)
    msg <- sprintf(
      "Detected %d row(s) with NA in: %s",
      n_bad, paste(cols_check, collapse = ", ")
    )
    if (na_action == "warn") {
      warning(msg)
    } else if (na_action == "drop") {
      df <- df[!na_rows, , drop = FALSE]
      message("Dropped ", n_bad, " row(s) with missing values")
    }
    # if na_action == "none": do nothing, but later operations may fail
  }


  # Both standardisation and relative fitness are computed from the rows that
  # are analysed (complete fitness and traits), within each group when
  # `group` is set. Individuals dropped from the gradient models for a missing
  # value would otherwise shift the trait mean/SD and the mean fitness and bias
  # every gradient (Lande & Arnold standardise within the sample under
  # selection). Rows whose group label is missing are kept and treated as their
  # own group; use na_action = "drop" to exclude them. The arithmetic is in
  # the internal utilities below, which bootstrap_selection() also runs on
  # every resample.
  cc <- stats::complete.cases(df[, c(fitness_col, trait_cols), drop = FALSE])
  grp_key <- if (!is.null(group)) addNA(factor(df[[group]]), ifany = TRUE) else rep(1L, nrow(df))
  groups <- .group_rows(if (!is.null(group)) as.integer(grp_key), nrow(df))
  group_label <- function(rows) {
    paste(unique(as.character(grp_key[rows])), collapse = ", ")
  }

  # Standardise traits, z = (x - mean(x)) / sd(x). A trait with no variance in
  # a group (constant, or a single observation) cannot be scaled there: it is
  # centred only, so those rows stay in the analysis with z = 0 rather than
  # becoming NaN, and the other groups are still scaled normally. A trait with
  # no variance in any group is left unscaled and named in a warning.
  if (standardize) {
    for (t in trait_cols) {
      s <- .scaling_within(df[[t]], cc, groups)
      if (!any(s$scalable[cc])) {
        warning("Zero-variance trait left unstandardised: ", t)
        next
      }
      if (any(!s$scalable[cc])) {
        warning(
          "Trait '", t, "' has no variance in group(s) ",
          group_label(cc & !s$scalable), "; centred only there"
        )
      }
      df[[t]] <- .standardise(df[[t]], s)
    }
  }

  # Add relative fitness, w_i = W_i / mean(W), with mean(W) over the analysed
  # rows. A group whose mean fitness is zero (nobody survived) has no relative
  # fitness: those rows are set to NA and named in a warning as they drop out
  # of the models.
  if (add_relative) {
    m <- .mean_fitness_within(df[[fitness_col]], cc, groups)
    if (!any(m$usable[cc])) {
      warning("Mean fitness is zero or non-finite; cannot compute relative fitness. Skipping.")
    } else {
      if (any(!m$usable[cc])) {
        warning(
          "Mean fitness is zero or non-finite in group(s) ", group_label(cc & !m$usable),
          "; relative fitness is NA there and those rows are excluded from the models"
        )
      }
      df[[name_relative]] <- .relativise(df[[fitness_col]], m)
    }
  }

  return(df)
}

#' @noRd
# internal utility: the rows of each group, from integer group codes; one
# group of every row without them
.group_rows <- function(codes, n) {
  if (is.null(codes)) list(seq_len(n)) else split(seq_len(n), codes)
}

#' @noRd
# internal utility: each row's statistic FUN of x over the analysed rows (cc)
# of its group, as stats::ave() gives it, written into a copy of the column.
# The ifelse() decides whether an integer column with nothing missing is
# averaged as integers, which can differ in the last bit from averaging it as
# doubles; and a column that is not numeric keeps its type, so its statistic
# is not finite.
.group_stat <- function(x, cc, groups, FUN) {
  x <- ifelse(cc, x, NA_real_)
  out <- x
  for (i in groups) {
    v <- x[i]
    out[i] <- FUN(v, na.rm = TRUE)
  }
  out
}

#' @noRd
# internal utility: each row's within-group mean and SD of x over the
# analysed rows, and `scalable`, FALSE for the rows of a group where x has no
# variance
.scaling_within <- function(x, cc, groups) {
  mu <- .group_stat(x, cc, groups, mean)
  sdev <- .group_stat(x, cc, groups, stats::sd)
  list(mu = mu, sdev = sdev, scalable = is.finite(sdev) & sdev > 0)
}

#' @noRd
# internal utility: z = (x - mean) / sd with the statistics from
# .scaling_within(); a group where x has no variance is only centred
.standardise <- function(x, s) {
  (x - s$mu) / ifelse(s$scalable, s$sdev, 1)
}

#' @noRd
# internal utility: each row's within-group mean fitness over the analysed
# rows, and `usable`, FALSE for the rows of a group whose mean is zero or not
# finite
.mean_fitness_within <- function(w, cc, groups) {
  mean_fit <- .group_stat(w, cc, groups, mean)
  list(mean = mean_fit, usable = is.finite(mean_fit) & mean_fit != 0)
}

#' @noRd
# internal utility: relative fitness w / mean(w) with the means from
# .mean_fitness_within(); NA in a group without a usable mean
.relativise <- function(w, m) {
  ifelse(m$usable, w / m$mean, NA_real_)
}
