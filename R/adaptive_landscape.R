# ======================================================
# adaptive_landscape.R
# Calculate Adaptive Landscape (Sewall Wright's concept)
#
# IMPORTANT CONCEPT:
# Adaptive Landscape = Mean fitness ~ Population mean phenotype
# Formula: Wbar ~ zbar1 + zbar2
#
# This is DIFFERENT from correlated fitness surface:
# - Correlated fitness: w ~ z (individual fitness)
# - Adaptive landscape: Wbar ~ zbar (population mean fitness)
#
# KEY PRINCIPLE:
# - Traits MUST already be standardized (mean = 0, SD = 1)
# - Use prepare_selection_data() with standardize = TRUE before calling
# - DO NOT standardize again within this function
#
# Requires:
# - A fitted fitness model (from GAM or TPS) that predicts individual fitness
# - Individual-level data to estimate within-population variance
# ======================================================

#' Calculate Adaptive Landscape
#'
#' Computes the adaptive landscape (mean fitness as a function of population mean phenotype)
#' using a fitted individual-level fitness model.
#'
#' @param data A data frame containing the original trait and fitness data.
#' @param fitness_model A fitted model object (e.g., GAM or Tps) predicting individual fitness.
#' @param trait_cols One or two trait column names. With one trait the model
#'   is usually the \code{$model} of \code{univariate_spline()} and the result
#'   is a curve; with two it is the \code{$model} of
#'   \code{correlated_fitness_surface()} and the result is a surface.
#' @param group_col Optional character string specifying a grouping variable.
#' @param population_variance Optional covariance matrix for the traits. Estimated from data if \code{NULL}.
#' @param simulation_n Integer specifying the number of individuals to simulate per grid point. Default is 1000.
#' @param grid_n Integer specifying the resolution of the population mean grid. Default is 50.
#' @param custom_range Optional list specifying custom ranges for the traits.
#' @param clamp Logical; if \code{TRUE} (the default) simulated fitness from a
#'   thin-plate model is held inside the range of the fitness type before
#'   averaging: 0 to 1 for survival, at least 0 for counts. Continuous fitness
#'   is left alone, and a GAM is unaffected because its link already respects
#'   the range. The run message says how many values were held.
#' @param support_warn Share of the population simulated at the optimum that
#'   may fall outside the data before the optimum is flagged as
#'   extrapolated. Default is 0.25.
#' @details For a single trait the grid also carries \code{.ind_fit}, the
#'   individual fitness function evaluated at each population mean, for
#'   drawing the two curves together (see \code{plot_adaptive_landscape()}).
#'
#'   The simulated populations spread beyond the data, especially towards the
#'   edge of the grid, and the fitness of those individuals is extrapolated.
#'   \code{.outside} in the grid is the
#'   share of each simulated population falling outside the convex hull of the
#'   observed trait pairs, or outside the observed range for one trait, and the
#'   result's \code{support} gives that share over the whole grid and at the
#'   optimum.
#' @return An object of class \code{"adaptive_landscape"}. Its
#'   \code{optimum_edge} is \code{TRUE} when the highest mean fitness lies on
#'   the edge of the grid, where the landscape may keep rising beyond it.
#' @examples
#' prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
#' surf <- correlated_fitness_surface(prep, "survival", c("total_length", "weight"), grid_n = 30)
#' land <- adaptive_landscape(prep, surf$model, c("total_length", "weight"),
#'                            simulation_n = 100, grid_n = 15)
#' land$optimum
#'
#' # one trait, from the spline fitness function
#' uni <- univariate_spline(prep, "survival", "total_length")
#' land1 <- adaptive_landscape(prep, uni$model, "total_length", simulation_n = 100, grid_n = 30)
#' plot_adaptive_landscape(land1, "total_length")
#' @export
adaptive_landscape <- function(
  data,
  fitness_model,
  trait_cols,
  group_col = NULL,
  population_variance = NULL,
  simulation_n = 1000,
  grid_n = 50,
  custom_range = NULL,
  clamp = TRUE,
  support_warn = 0.25
) {
    # Input validation
    stopifnot(length(trait_cols) %in% c(1L, 2L))
    stopifnot(inherits(fitness_model, "gam") || inherits(fitness_model, "Tps"))

    # the traits come in standardised; warn when one does not look it
    for (t in trait_cols) {
        if (t %in% names(data)) {
            z_mean <- mean(data[[t]], na.rm = TRUE)
            z_sd <- sd(data[[t]], na.rm = TRUE)
            if (abs(z_mean) > 0.1 || abs(z_sd - 1) > 0.1) {
                warning(
                    "Trait '", t, "' does not look standardised ",
                    "(mean ", round(z_mean, 3), ", SD ", round(z_sd, 3), "), so it is used in its own units; ",
                    "prepare_selection_data() standardises it"
                )
            }
        }
    }


    # Range of the population mean grid, one per trait
    ranges <- lapply(trait_cols, function(t) {
        if (!is.null(custom_range) && !is.null(custom_range[[t]])) {
            return(custom_range[[t]])
        }
        # Expand range slightly to include possible evolutionary space
        r <- range(data[[t]], na.rm = TRUE)
        r + c(-0.2, 0.2) * diff(r)
    })
    names(ranges) <- trait_cols

    message("Population mean grid ranges:\n", paste0(
        "  ", trait_cols, ": ",
        vapply(trait_cols, function(t) paste(round(ranges[[t]], 2), collapse = " to "), character(1)),
        collapse = "\n"
    ))

    # Create grid of population mean phenotypes
    population_grid <- expand.grid(
        lapply(ranges, function(r) seq(r[1], r[2], length.out = grid_n)),
        KEEP.OUT.ATTRS = FALSE
    )
    names(population_grid) <- trait_cols

    # Estimate within-population variance if not provided
    if (is.null(population_variance)) {
        # Use all data to estimate phenotypic variance
        trait_all <- data[, trait_cols, drop = FALSE]
        cc <- complete.cases(trait_all)
        trait_data <- trait_all[cc, , drop = FALSE]

        if (nrow(trait_data) > 1) {
            population_variance <- var(trait_data, na.rm = TRUE)
        } else {
            stop("Not enough data to estimate population variance")
        }

        # If group_col provided, use average within-group variance
        if (!is.null(group_col) && group_col %in% names(data)) {
            group_variances <- by(
                trait_data,
                data[[group_col]][cc],
                function(x) if (nrow(x) > 1) var(x) else NULL
            )
            group_variances <- group_variances[!sapply(group_variances, is.null)]
            if (length(group_variances) > 0) {
                population_variance <- Reduce(`+`, group_variances) / length(group_variances)
                message("Using average within-group variance")
            }
        }

        .msg_table("Estimated within-population variance-covariance:", round(population_variance, 4))
    }

    # Ensure variance matrix is positive definite
    if (inherits(try(chol(population_variance), silent = TRUE), "try-error")) {
        message("Making the variance matrix positive definite")
        if (requireNamespace("Matrix", quietly = TRUE)) {
            population_variance <- as.matrix(Matrix::nearPD(population_variance)$mat)
        } else {
            # Simple correction: add small constant to diagonal
            population_variance <- population_variance + diag(1e-6, nrow(population_variance))
        }
    }

    mean_fitness <- numeric(nrow(population_grid))

    # only the thin-plate model can predict outside the fitness range
    do_clamp <- isTRUE(clamp) && inherits(fitness_model, "Tps")
    ftype <- if (do_clamp) suppressWarnings(detect_family(as.numeric(fitness_model$y)))$type else NULL
    clipped <- 0L

    # how much of each simulated population leaves the data: outside the convex
    # hull of the observed trait pairs, or outside the observed range for one trait
    obs <- data[stats::complete.cases(data[, trait_cols, drop = FALSE]), trait_cols, drop = FALSE]
    hull <- if (length(trait_cols) == 2L) .data_hull(obs[[1]], obs[[2]]) else NULL
    obs_range <- range(obs[[1]])
    outside <- numeric(nrow(population_grid))

    # A GAM fitted with a group term (correlated_fitness_surface(..., group =))
    # needs that variable to predict from. Hold every non-trait predictor at a
    # reference level taken from the data the model was fitted to, where the
    # group is a factor, so a numeric group such as year is held at the same
    # level as on the surface (as peak_difference() does).
    extra_vars <- list()
    if (inherits(fitness_model, "gam")) {
        mf_names <- names(fitness_model$model)
        model_vars <- setdiff(mf_names[-1], trait_cols) # first column is the response
        model_vars <- model_vars[!grepl("^\\(", model_vars)] # skip (weights), (offset)
        for (v in model_vars) {
            if (!v %in% names(data)) {
                stop("The fitness model uses '", v, "', which is not in `data`")
            }
            extra_vars[[v]] <- .reference_group(fitness_model$model[[v]])
            message("Holding '", v, "' at reference level '", extra_vars[[v]], "'")
        }
    }

    # Progress indicator
    message("Calculating mean fitness for ", nrow(population_grid), " grid points")

    for (i in seq_len(nrow(population_grid))) {
        pop_mean <- population_grid[i, ]

        # Simulate individuals around this population mean
        simulated <- MASS::mvrnorm(
            n = simulation_n,
            mu = as.numeric(pop_mean),
            Sigma = population_variance
        )

        colnames(simulated) <- trait_cols
        outside[i] <- if (is.null(hull)) {
            mean(simulated[, 1] < obs_range[1] | simulated[, 1] > obs_range[2])
        } else {
            mean(!mgcv::in.out(hull, unname(simulated[, 1:2, drop = FALSE])))
        }
        sim_df <- as.data.frame(simulated)
        for (v in names(extra_vars)) sim_df[[v]] <- extra_vars[[v]]

        # Predict individual fitness
        if (inherits(fitness_model, "gam")) {
            ind_fitness <- predict(fitness_model, newdata = sim_df, type = "response")
        } else if (inherits(fitness_model, "Tps")) {
            ind_fitness <- predict(fitness_model, as.matrix(sim_df))
        }
        if (do_clamp) {
            held <- .clamp_fitness(ind_fitness, ftype)
            clipped <- clipped + sum(held != ind_fitness, na.rm = TRUE)
            ind_fitness <- held
        }

        # Mean population fitness
        mean_fitness[i] <- mean(ind_fitness, na.rm = TRUE)
    }


    if (do_clamp && clipped > 0) {
        message("Held ", clipped, " of ", simulation_n * nrow(population_grid),
                " simulated fitness values inside the range of ", ftype, " fitness")
    }

    # Add mean fitness to grid
    population_grid$.mean_fit <- mean_fitness
    population_grid$.outside <- outside

    # With one trait keep the individual fitness function at the same means,
    # so the landscape can be drawn against it.
    if (length(trait_cols) == 1L) {
        ind_df <- population_grid[, trait_cols, drop = FALSE]
        for (v in names(extra_vars)) ind_df[[v]] <- extra_vars[[v]]
        population_grid$.ind_fit <- if (inherits(fitness_model, "gam")) {
            as.numeric(predict(fitness_model, newdata = ind_df, type = "response"))
        } else {
            as.numeric(predict(fitness_model, as.matrix(ind_df)))
        }
        if (do_clamp) population_grid$.ind_fit <- .clamp_fitness(population_grid$.ind_fit, ftype)
    }

    # Find optimum (maximum mean fitness)
    optimum <- population_grid[which.max(mean_fitness), ]
    .msg_table("Optimal population mean phenotype:", optimum[, trait_cols, drop = FALSE])
    message("Mean fitness at optimum: ", round(optimum$.mean_fit, 4))
    # a highest point on the edge of the grid may not be an optimum
    optimum_edge <- any(vapply(trait_cols, function(t) {
        v <- population_grid[[t]]
        optimum[[t]] <= min(v) + 1e-9 || optimum[[t]] >= max(v) - 1e-9
    }, logical(1)))
    if (optimum_edge) {
        message("The highest mean fitness is on the edge of the grid; the landscape may keep rising beyond it")
    }

    support <- list(outside = mean(outside), at_optimum = optimum$.outside, warn = support_warn)
    message(sprintf("%.0f%% of simulated individuals fell outside the data (%.0f%% at the optimum)",
                    100 * support$outside, 100 * support$at_optimum))
    if (support$at_optimum > support_warn) {
        message("The optimum rests on extrapolation: more than ", 100 * support_warn,
                "% of the population simulated there lies outside the data")
    }

    # Calculate actual population means if group_col provided
    actual_means <- NULL
    if (!is.null(group_col) && group_col %in% names(data)) {
        actual_means <- aggregate(
            data[, trait_cols, drop = FALSE],
            by = list(group = data[[group_col]]),
            FUN = mean,
            na.rm = TRUE
        )
        names(actual_means)[1] <- group_col
        .msg_table("Actual population means:", actual_means)
    }

    result <- list(
        grid = population_grid,
        trait_cols = trait_cols,
        population_variance = population_variance,
        simulation_n = simulation_n,
        clamp = do_clamp,
        clipped = clipped,
        support = support,
        optimum = optimum,
        optimum_edge = optimum_edge,
        actual_population_means = actual_means,
        fitness_model_class = class(fitness_model)[1],
        data_summary = list(
            n_individuals = nrow(data),
            trait_ranges = stats::setNames(
                lapply(trait_cols, function(t) range(data[[t]], na.rm = TRUE)),
                paste0("x", seq_along(trait_cols))
            ),
            traits_standardized = all(
                sapply(trait_cols, function(t) {
                    abs(mean(data[[t]], na.rm = TRUE)) < 0.1 &&
                        abs(sd(data[[t]], na.rm = TRUE) - 1) < 0.1
                })
            )
        ),
        surface_type = "adaptive_landscape",
        note = "Adaptive landscape: mean fitness as function of population mean phenotype"
    )

    class(result) <- "adaptive_landscape"
    return(result)
}


#' Print method for adaptive landscape
#'
#' @param x An object of class \code{"adaptive_landscape"}.
#' @param ... Additional arguments passed to \code{print}.
#' @return The input object \code{x}, invisibly.
#' @export
print.adaptive_landscape <- function(x, ...) {
    cat("\nAdaptive Landscape Object\n")
    cat("========================\n")
    cat(if (length(x$trait_cols) == 1L) "Trait:" else "Traits:",
        paste(x$trait_cols, collapse = " and "), "\n")
    cat("Grid size:", nrow(x$grid), "population means\n")
    cat("Simulations per point:", x$simulation_n, "individuals\n")
    cat("Fitness model:", x$fitness_model_class, "\n")

    # Check if traits were standardized
    if (!is.null(x$data_summary$traits_standardized)) {
        if (x$data_summary$traits_standardized) {
            cat("Traits: standardised (mean ~ 0, SD ~ 1)\n")
        } else {
            cat("Traits: may not be standardised\n")
        }
    }

    cat("\nOptimal population mean phenotype:\n")
    print(x$optimum[, x$trait_cols, drop = FALSE])
    cat("Mean fitness at optimum:", round(x$optimum$.mean_fit, 4), "\n")
    if (isTRUE(x$optimum_edge)) cat("It lies on the edge of the grid; the landscape may keep rising beyond it\n")
    if (!is.null(x$support)) {
        cat(sprintf("Simulated individuals outside the data: %.0f%% over the grid, %.0f%% at the optimum\n",
                    100 * x$support$outside, 100 * x$support$at_optimum))
        if (x$support$at_optimum > x$support$warn) cat("The optimum rests on extrapolation beyond the data\n")
    }

    if (!is.null(x$actual_population_means)) {
        cat("\nActual population means:\n")
        print(x$actual_population_means)
    }
    cat("\n", x$note, "\n")
    invisible(x)
}

# hold predictions inside the range of the fitness type: 0 to 1 for survival,
# at least 0 for counts; continuous fitness is left alone
.clamp_fitness <- function(x, type) {
    if (identical(type, "binary")) pmin(pmax(x, 0), 1)
    else if (identical(type, "count")) pmax(x, 0)
    else x
}
