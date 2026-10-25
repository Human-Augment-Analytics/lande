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
#'   thin-plate model (with \code{"convolve"}, its fitness on the lattice) is
#'   held inside the range of the fitness type before
#'   averaging: 0 to 1 for survival, at least 0 for counts. Continuous fitness
#'   is left alone, and a GAM is unaffected because its link already respects
#'   the range. The run message says how many values were held.
#' @param support_warn Share of the population simulated at the optimum that
#'   may fall outside the data before the optimum is flagged as
#'   extrapolated. Default is 0.25.
#' @param method How the mean fitness of each population is worked out:
#'   \code{"simulate"} (the default), \code{"lattice"} or \code{"convolve"}.
#'   See Details.
#' @param lattice_n Nodes per trait of the lattice the \code{"lattice"}
#'   method predicts the fitness surface on, for two traits (at least 8). One
#'   trait uses \code{max(10 * lattice_n, 1000)} nodes. Default is 200.
#' @details All three methods estimate the same quantity, the mean fitness of
#'   a population whose traits are normally distributed around each grid
#'   point with the within-population (co)variance:
#'   \itemize{
#'     \item \code{"simulate"} draws \code{simulation_n} individuals at every
#'       grid point, predicts the fitness of each one from the model and
#'       averages them (Monte Carlo). It is the slowest, and the answer
#'       moves a little with the random seed.
#'     \item \code{"lattice"} makes the same random draws, in the same order,
#'       but predicts the fitness surface once on a fine lattice and
#'       interpolates it for every individual (cubic interpolation along
#'       each trait). For a seed it gives nearly the
#'       same numbers as \code{"simulate"}, several times faster, and leaves
#'       the random number stream where \code{"simulate"} would. A sample of
#'       the draws is also predicted exactly to check the interpolation; the
#'       lattice is made finer if the error is more than 0.1\% of the range
#'       of fitness over the grid of population means, and the result's
#'       \code{lattice_check} records it. When the lattice would need as many
#'       predictions as the draws, or a finer one would cost more than
#'       predicting the remaining draws, the draws are predicted directly as
#'       in \code{"simulate"}.
#'     \item \code{"convolve"} draws nothing. It evaluates the fitness surface
#'       on a lattice that contains every grid point and averages it with the
#'       weights of the normal distribution (a convolution), so there is no
#'       simulation noise, the landscape is smooth and the result does not
#'       depend on the seed. \code{simulation_n} is not used. The lattice is
#'       at most a quarter of each trait's SD apart, and finer when the
#'       traits are strongly correlated, and the weights reach 5 SD around
#'       each grid point (8 SD for a GAM with a log link, whose fitness can
#'       climb steeply). Its cost depends on the grid and the variance, not
#'       on \code{simulation_n}: when the distribution is very narrow next to
#'       the grid step, or the traits are almost perfectly correlated, it
#'       would need more predictions than simulating \code{simulation_n}
#'       individuals at every grid point, and the \code{"lattice"} method is
#'       used instead (with a message).
#'   }
#'   The \code{"lattice"} method interpolates a GAM on the scale of its link
#'   function and transforms the result back, so fitness stays inside the
#'   range the link allows.
#'
#'   For a single trait the grid also carries \code{.ind_fit}, the
#'   individual fitness function evaluated at each population mean, for
#'   drawing the two curves together (see \code{plot_adaptive_landscape()}).
#'
#'   The simulated populations spread beyond the data, especially towards the
#'   edge of the grid, and the fitness of those individuals is extrapolated.
#'   \code{.outside} in the grid is the
#'   share of each simulated population falling outside the convex hull of the
#'   observed trait pairs, or outside the observed range for one trait, and the
#'   result's \code{support} gives that share over the whole grid and at the
#'   optimum. With \code{"convolve"} it is the share of each population's
#'   normal distribution lying outside the data, taking the part of each
#'   lattice cell that lies outside.
#' @return An object of class \code{"adaptive_landscape"}. Its
#'   \code{optimum_edge} is \code{TRUE} when the highest mean fitness lies on
#'   the edge of the grid, where the landscape may keep rising beyond it.
#'   \code{method} is the method used (which can differ from the one asked
#'   for when it falls back), \code{lattice} describes the lattice
#'   of the \code{"lattice"} and \code{"convolve"} methods (nodes and spacing
#'   per trait, number of predictions) and \code{lattice_check} the
#'   interpolation check of the \code{"lattice"} method: draws predicted
#'   exactly, the largest error among them, the tolerance, whether it
#'   passed, the final nodes per trait, the number of refinements and the
#'   number of draws predicted directly. \code{clipped} counts the fitness
#'   values held by \code{clamp}: simulated individuals for
#'   \code{"simulate"} and \code{"lattice"}, lattice nodes for
#'   \code{"convolve"}.
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
#'
#' # the same landscape without simulation noise
#' land_c <- adaptive_landscape(prep, surf$model, c("total_length", "weight"),
#'                              grid_n = 15, method = "convolve")
#' land_c$optimum
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
  support_warn = 0.25,
  method = c("simulate", "lattice", "convolve"),
  lattice_n = 200
) {
    # Input validation
    stopifnot(length(trait_cols) %in% c(1L, 2L))
    stopifnot(inherits(fitness_model, "gam") || inherits(fitness_model, "Tps"))
    method <- match.arg(method)
    # the convolution falls back to the lattice when it would cost too much
    if (method != "simulate") {
        stopifnot(is.numeric(lattice_n), length(lattice_n) == 1L, is.finite(lattice_n), lattice_n >= 8)
        lattice_n <- as.integer(round(lattice_n))
    }
    simulation_n_given <- !missing(simulation_n)
    # a model read back from disk predicts through its own package, so load it
    # and its predict() method before anything is predicted
    invisible(requireNamespace(if (inherits(fitness_model, "gam")) "mgcv" else "fields", quietly = TRUE))

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

    lattice <- NULL
    lattice_check <- NULL
    # the convolution's lattice grows as the distribution narrows in some
    # direction; past the predictions the simulation would make, interpolating
    # simulated individuals is the better deal
    if (method == "convolve") {
        # a log link lets fitness climb steeply, so its weights reach further
        reach <- if (inherits(fitness_model, "gam") && identical(fitness_model$family$link, "log")) 8 else 5
        plan <- .convolve_plan(population_grid, trait_cols, population_variance, reach)
        budget <- max(simulation_n * nrow(population_grid), 5e4)
        if (plan$predictions > budget) {
            count <- function(x) {
                if (!is.finite(x)) "far too many"
                else if (x < 1e9) format(round(x), big.mark = ",", scientific = FALSE)
                else format(signif(x, 2))
            }
            message("Convolution would need ", count(plan$predictions), " fitness predictions here, more than the ",
                    count(budget), " allowed (as many as the simulated individuals, and at least 50,000), as the ",
                    "populations' distribution is very narrow next to the grid step or the traits are almost ",
                    "perfectly correlated; using the lattice method instead")
            method <- "lattice"
        } else if (simulation_n_given) {
            message("The convolve method draws no individuals, so simulation_n is not used")
        }
    }
    # a lattice is only worth it when it needs fewer predictions than the draws
    if (method == "lattice") {
        n_draws <- simulation_n * nrow(population_grid)
        n_nodes <- .lattice_axis_n(lattice_n, length(trait_cols))^length(trait_cols)
        if (n_nodes >= n_draws) {
            message("The lattice would need ", format(n_nodes, scientific = FALSE), " predictions, no fewer than the ",
                    format(n_draws, scientific = FALSE),
                    " simulated individuals, so each one is predicted directly (method \"simulate\")")
            method <- "simulate"
        }
    }

    if (method == "lattice") {
        lat <- .landscape_lattice(fitness_model, population_grid, trait_cols, population_variance,
                                  simulation_n, lattice_n, extra_vars, do_clamp, ftype, hull, obs_range)
        mean_fitness <- lat$mean_fit
        outside <- lat$outside
        clipped <- lat$clipped
        lattice <- lat$lattice
        lattice_check <- lat$check
    } else if (method == "convolve") {
        conv <- .landscape_convolve(fitness_model, plan, trait_cols, population_variance,
                                    extra_vars, do_clamp, ftype, hull, obs_range)
        mean_fitness <- conv$mean_fit
        outside <- conv$outside
        clipped <- conv$clipped
        lattice <- conv$lattice
    } else {
        # Monte Carlo: draw a population at every grid point and predict each individual
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
    }


    if (do_clamp && clipped > 0) {
        if (method == "convolve") {
            message("Held ", clipped, " of ", lattice$predictions,
                    " fitness values on the lattice inside the range of ", ftype, " fitness")
        } else {
            message("Held ", clipped, " of ", simulation_n * nrow(population_grid),
                    " simulated fitness values inside the range of ", ftype, " fitness")
        }
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
    if (method == "convolve") {
        message(sprintf("%.0f%% of the populations' distribution lay outside the data (%.0f%% at the optimum)",
                        100 * support$outside, 100 * support$at_optimum))
    } else {
        message(sprintf("%.0f%% of simulated individuals fell outside the data (%.0f%% at the optimum)",
                        100 * support$outside, 100 * support$at_optimum))
    }
    if (support$at_optimum > support_warn) {
        message("The optimum rests on extrapolation: more than ", 100 * support_warn,
                if (method == "convolve") "% of the population there lies outside the data"
                else "% of the population simulated there lies outside the data")
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
        # the convolution draws no individuals
        simulation_n = if (method == "convolve") NA_integer_ else simulation_n,
        method = method,
        lattice = lattice,
        lattice_check = lattice_check,
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
    # objects made before the methods existed were simulated
    method <- if (is.null(x$method)) "simulate" else x$method
    if (method == "convolve") {
        cat("Method: convolve (the population's normal distribution averaged over a",
            paste(x$lattice$n, collapse = " x "), "lattice, no simulation)\n")
    } else {
        cat("Method:", if (method == "lattice") "lattice (simulated individuals, fitness interpolated from a lattice)"
            else "simulate (simulated individuals)", "\n")
        cat("Simulations per point:", x$simulation_n, "individuals\n")
    }
    if (method == "lattice" && !is.null(x$lattice_check)) {
        chk <- x$lattice_check
        if (isTRUE(chk$direct > 0)) {
            cat(chk$direct, "draws predicted directly, the lattice not being fine enough\n")
        }
        if (!isTRUE(chk$direct >= nrow(x$grid) * x$simulation_n)) {
            cat(sprintf("Lattice check: largest error %.2g at %d draws predicted exactly (tolerance %.2g)\n",
                        chk$max_error, chk$n_checked, chk$tol))
        }
    }
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
        cat(sprintf("%s outside the data: %.0f%% over the grid, %.0f%% at the optimum\n",
                    if (method == "convolve") "Population" else "Simulated individuals",
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

# nodes per trait of the interpolation lattice: lattice_n for two traits, and
# a cheap fine line of at least 1000 for one
.lattice_axis_n <- function(lattice_n, p) {
    if (p == 1L) max(10L * lattice_n, 1000L) else lattice_n
}

# predict the fitness model at the rows of `pts` (one column per trait),
# holding any other predictor at its reference level, on the scale of the
# GAM's linear predictor when link = TRUE. Done in blocks so a fine lattice
# does not build one huge matrix.
.predict_fitness <- function(fitness_model, pts, trait_cols, extra_vars, link = FALSE, block = 50000L) {
    pts <- matrix(pts, ncol = length(trait_cols), dimnames = list(NULL, trait_cols))
    out <- numeric(nrow(pts))
    rows_all <- seq_len(nrow(pts))
    for (rows in split(rows_all, ceiling(rows_all / block))) {
        df <- as.data.frame(pts[rows, , drop = FALSE])
        for (v in names(extra_vars)) df[[v]] <- extra_vars[[v]]
        out[rows] <- if (inherits(fitness_model, "gam")) {
            as.numeric(predict(fitness_model, newdata = df, type = if (link) "link" else "response"))
        } else {
            as.numeric(predict(fitness_model, as.matrix(df)))
        }
    }
    out
}

# interpolation of the values on a regular lattice at the rows of X: cubic
# (Catmull-Rom) along each trait, which follows a sharply bending surface far
# more closely than linear interpolation on the same lattice. The lattice
# values run with the first trait fastest, as expand.grid() lays them out. NA
# where a row falls off the lattice or in its outermost cells, where the
# 4 nodes around it along a trait are not all on the lattice.
.interp_lattice <- function(lat, X) {
    nx <- lat$n[1]
    fx <- (X[, 1] - lat$lo[1]) / lat$h[1]
    out <- rep(NA_real_, nrow(X))
    if (length(lat$n) == 1L) {
        ix <- floor(fx)
        ok <- !is.na(ix) & ix >= 1 & ix <= nx - 3
        ix <- ix[ok]
        wx <- .catmull_rom(fx[ok] - ix)
        acc <- 0
        for (a in 1:4) acc <- acc + wx[[a]] * lat$z[ix + a - 1]
        out[ok] <- acc
        return(out)
    }
    ny <- lat$n[2]
    fy <- (X[, 2] - lat$lo[2]) / lat$h[2]
    ix <- floor(fx)
    iy <- floor(fy)
    # the 4 x 4 nodes around a row must all be on the lattice
    ok <- !is.na(ix) & !is.na(iy) & ix >= 1 & ix <= nx - 3 & iy >= 1 & iy <= ny - 3
    ix <- ix[ok]
    iy <- iy[ok]
    wx <- .catmull_rom(fx[ok] - ix)
    wy <- .catmull_rom(fy[ok] - iy)
    k <- ix + 1 + iy * nx
    acc <- 0
    for (b in 1:4) {
        along <- 0
        for (a in 1:4) along <- along + wx[[a]] * lat$z[k + (a - 2) + (b - 2) * nx]
        acc <- acc + wy[[b]] * along
    }
    out[ok] <- acc
    out
}

# Catmull-Rom weights of the nodes at -1, 0, 1 and 2 for a point a fraction t
# of the way from node 0 to node 1; they sum to 1 and give the node values at
# the nodes themselves
.catmull_rom <- function(t) {
    t2 <- t * t
    t3 <- t2 * t
    list((-t3 + 2 * t2 - t) / 2, (3 * t3 - 5 * t2 + 2) / 2, (-3 * t3 + 4 * t2 + t) / 2, (t3 - t2) / 2)
}

# method = "lattice": the same random draws as the simulate loop, in the same
# order, but the fitness of each individual is interpolated from the surface
# predicted once on a lattice over the grid and 5 SD beyond it. A fixed
# stride of the draws is predicted exactly to check the interpolation (no
# random numbers). When the check fails the lattice is made finer, if one
# fine enough to pass should cost fewer predictions than the draws left;
# otherwise the draws left are predicted directly.
.landscape_lattice <- function(fitness_model, population_grid, trait_cols, population_variance,
                               simulation_n, lattice_n, extra_vars, do_clamp, ftype, hull, obs_range) {
    p <- length(trait_cols)
    n_grid <- nrow(population_grid)
    n_draws <- simulation_n * n_grid
    sds <- sqrt(diag(as.matrix(population_variance)))
    lo <- vapply(trait_cols, function(t) min(population_grid[[t]]), numeric(1)) - 5 * sds
    hi <- vapply(trait_cols, function(t) max(population_grid[[t]]), numeric(1)) + 5 * sds

    # a GAM is interpolated on the scale of its linear predictor and
    # transformed back, so the fitness keeps to the range its link allows
    linkinv <- if (inherits(fitness_model, "gam")) fitness_model$family$linkinv else NULL
    on_link <- is.function(linkinv)
    to_response <- function(x) if (on_link) linkinv(x) else x
    hold <- function(x) if (do_clamp) .clamp_fitness(x, ftype) else x

    # the check allows a thousandth of the range of fitness over the grid of
    # population means, which is about the range of the landscape itself.
    # Fitness out in the lattice's margin can be far larger, under a log
    # link, and would loosen the check.
    spread <- function(v) {
        v <- v[is.finite(v)]
        if (length(v) >= 2L) diff(range(v)) else 0
    }
    tol_grid <- 1e-3 * spread(hold(.predict_fitness(fitness_model, as.matrix(population_grid[, trait_cols]),
                                                    trait_cols, extra_vars)))

    make_lattice <- function(n_axis) {
        n <- rep(n_axis, p)
        nodes <- as.matrix(expand.grid(lapply(seq_len(p), function(j) seq(lo[j], hi[j], length.out = n[j])),
                                       KEEP.OUT.ATTRS = FALSE))
        z <- .predict_fitness(fitness_model, nodes, trait_cols, extra_vars, link = on_link)
        # a flat grid (or a single population) has no range to go by: take the lattice's
        tol <- tol_grid
        if (!(tol > 0)) {
            resp <- hold(to_response(z))
            tol <- 1e-3 * spread(resp)
            if (!(tol > 0) || !is.finite(tol)) tol <- 1e-8 * max(1, abs(resp[is.finite(resp)]))
        }
        list(n = n, lo = lo, h = (hi - lo) / (n - 1), z = z, tol = tol)
    }
    lat_size <- function(lat) paste(lat$n, collapse = " x ")

    # nodes per trait at most: 800 for two traits, a long line for one
    cap <- max(if (p == 1L) 1e5 else 800, .lattice_axis_n(lattice_n, p))
    # how many times the error falls when the spacing halves: about 8 for
    # cubic interpolation once the lattice is fine; a refinement that bought
    # less lowers it, though not below 4, as a coarse lattice can gain little
    # from one step and much from the next
    rate <- 8
    # the predictions of the finer lattices needed to bring the error down
    # by `ratio`, up to the cap
    refine_cost <- function(n_axis, ratio) {
        steps <- max(1, ceiling(log(ratio) / log(rate)))
        cost <- 0
        while (steps > 0 && n_axis < cap) {
            n_axis <- min(2 * n_axis, cap)
            cost <- cost + n_axis^p
            steps <- steps - 1
        }
        cost
    }

    lat <- make_lattice(.lattice_axis_n(lattice_n, p))
    predictions <- length(lat$z)
    refinements <- 0L
    exact_all <- FALSE
    n_exact <- 0L
    n_direct <- 0L
    max_error <- 0
    check_at <- seq(1L, n_draws, by = max(1L, ceiling(n_draws / 2000)))
    mean_fit <- outside <- numeric(n_grid)
    clipped <- 0L

    # whole grid points at a time, holding no more than about 250,000 draws
    per_chunk <- max(1L, floor(2.5e5 / simulation_n))
    for (idx in split(seq_len(n_grid), ceiling(seq_len(n_grid) / per_chunk))) {
        X <- do.call(rbind, lapply(idx, function(i) {
            sim <- MASS::mvrnorm(n = simulation_n, mu = as.numeric(population_grid[i, ]), Sigma = population_variance)
            # a single draw comes back as a vector
            if (is.matrix(sim)) sim else matrix(sim, ncol = p)
        }))
        off_data <- if (is.null(hull)) {
            X[, 1] < obs_range[1] | X[, 1] > obs_range[2]
        } else {
            !mgcv::in.out(hull, X[, 1:2, drop = FALSE])
        }
        outside[idx] <- colMeans(matrix(off_data, simulation_n))

        # the draws of this chunk that are in the check
        first <- (idx[1] - 1L) * simulation_n
        chk <- check_at[check_at > first & check_at <= first + nrow(X)] - first
        exact_chk <- if (exact_all) NULL else {
            hold(.predict_fitness(fitness_model, X[chk, , drop = FALSE], trait_cols, extra_vars))
        }

        prev_err <- NA_real_
        repeat {
            if (exact_all) {
                w <- .predict_fitness(fitness_model, X, trait_cols, extra_vars)
                n_direct <- n_direct + nrow(X)
                err <- 0
                break
            }
            z <- .interp_lattice(lat, X)
            # the rare draw beyond the lattice is predicted exactly
            off <- which(is.na(z))
            if (length(off)) {
                z[off] <- .predict_fitness(fitness_model, X[off, , drop = FALSE], trait_cols, extra_vars,
                                           link = on_link)
            }
            w <- to_response(z)
            err <- if (length(chk)) max(abs(hold(w[chk]) - exact_chk), na.rm = TRUE) else 0
            # what the last refinement actually bought
            if (!is.na(prev_err)) rate <- min(8, max(4, prev_err / err))
            # good enough, or as fine as the cap allows: keep it and report the check
            if (err <= lat$tol || lat$n[1] >= cap) {
                n_exact <- n_exact + length(union(chk, off))
                break
            }
            # a finer lattice pays only when it should pass for fewer
            # predictions than predicting the draws left directly
            if (refine_cost(lat$n[1], err / lat$tol) < n_draws - first) {
                message(sprintf("Lattice of %s nodes is off by %.2g (tolerance %.2g); trying a finer one",
                                lat_size(lat), err, lat$tol))
                prev_err <- err
                lat <- make_lattice(min(2 * lat$n[1], cap))
                predictions <- predictions + length(lat$z)
                refinements <- refinements + 1L
            } else {
                message(sprintf(paste0("Lattice of %s nodes is off by %.2g (tolerance %.2g) and one fine enough ",
                                       "would need more predictions than the draws left, so they are predicted directly"),
                                lat_size(lat), err, lat$tol))
                exact_all <- TRUE
            }
        }
        max_error <- max(max_error, err)

        if (do_clamp) {
            held <- .clamp_fitness(w, ftype)
            clipped <- clipped + sum(held != w, na.rm = TRUE)
            w <- held
        }
        mean_fit[idx] <- colMeans(matrix(w, simulation_n), na.rm = TRUE)
    }

    # max_error is that of the values used: nil for draws predicted directly
    check <- list(n_checked = length(check_at), max_error = max_error, tol = lat$tol,
                  passed = max_error <= lat$tol, lattice_n = lat$n[1], refinements = refinements,
                  direct = n_direct)
    if (n_direct == n_draws) {
        message(sprintf("Lattice of %s nodes was not fine enough, so all %d draws were predicted directly",
                        lat_size(lat), n_draws))
    } else if (n_direct > 0) {
        message(sprintf(paste0("Lattice of %s nodes: %d of the %d draws were predicted directly; the rest are off ",
                               "by at most %.2g at the draws predicted exactly, %s the tolerance %.2g"),
                        lat_size(lat), n_direct, n_draws, max_error,
                        if (check$passed) "within" else "above", lat$tol))
    } else {
        message(sprintf("Lattice of %s nodes: largest error %.2g at %d draws predicted exactly, %s the tolerance %.2g",
                        lat_size(lat), max_error, length(check_at),
                        if (check$passed) "within" else "above", lat$tol))
    }
    list(
        mean_fit = mean_fit,
        outside = outside,
        clipped = clipped,
        check = check,
        lattice = list(
            n = stats::setNames(lat$n, trait_cols),
            spacing = stats::setNames(lat$h, trait_cols),
            range = stats::setNames(lapply(seq_len(p), function(j) unname(c(lo[j], hi[j]))), trait_cols),
            predictions = predictions,
            # draws predicted exactly: the check, any beyond the lattice and any predicted directly
            exact = n_exact + n_direct
        )
    )
}

# method = "convolve", the plan: spacing, reach and number of nodes along each
# trait, worked out before anything is predicted so a lattice too large to be
# worth it is never built. The spacing is a whole fraction of the grid step,
# so every grid point is a node, at most a quarter of the trait's SD, and at
# most its SD with the other trait held fixed divided by sqrt(2). Strongly
# correlated traits squeeze the distribution into a narrow ridge running
# across the lattice, and the nodes must lie close enough to follow it. By
# Poisson summation the sum over the nodes is off from the integral by terms
# of size exp(-2 pi^2 k'H^-1 S H^-1 k) for whole k other than 0 (H the
# spacings), and each such k'H^-1 S H^-1 k is at least the smallest
# (conditional SD / spacing)^2, so this spacing keeps those terms below
# exp(-4 pi^2) however strong the correlation. The weights reach `reach` SDs
# along each trait.
.convolve_plan <- function(population_grid, trait_cols, population_variance, reach) {
    p <- length(trait_cols)
    S <- as.matrix(population_variance)
    sds <- sqrt(diag(S))
    cond_sds <- if (p == 2L) {
        inv <- tryCatch(solve(S), error = function(e) NULL)
        if (is.null(inv)) c(0, 0) else 1 / sqrt(diag(inv))
    } else {
        sds
    }
    h_max <- pmin(sds / 4, cond_sds / sqrt(2))
    n_g <- as.integer(round(nrow(population_grid)^(1 / p)))
    axes <- lapply(seq_len(p), function(j) {
        if (!is.finite(h_max[j]) || h_max[j] <= 0) return(list(n_need = Inf))
        # the grid's own values along this trait (the first trait runs fastest)
        g <- population_grid[[trait_cols[j]]][(seq_len(n_g) - 1L) * n_g^(j - 1L) + 1L]
        step <- if (n_g > 1L) (g[n_g] - g[1]) / (n_g - 1) else 0
        # r lattice steps to a grid step
        r <- if (step != 0) ceiling(abs(step) / h_max[j]) else 0
        h <- if (step != 0) step / r else h_max[j]
        m <- ceiling(reach * sds[j] / abs(h))
        # the nodes some grid point's weights reach: 2m + 1 around each grid
        # point, joining up when those overlap
        n_need <- if (r == 0) 2 * m + 1 else if (r <= 2 * m + 1) (n_g - 1) * r + 2 * m + 1 else n_g * (2 * m + 1)
        list(g1 = g[1], h = h, r = r, m = m, n_grid = n_g, n_need = n_need)
    })
    list(axes = axes, reach = reach,
         predictions = prod(vapply(axes, function(a) as.numeric(a$n_need), numeric(1))))
}

# method = "convolve": the mean fitness at each grid point as the fitness
# surface averaged with the weights of the population's normal distribution,
# a discrete convolution on a lattice that passes through every grid point and
# reaches 5 SD (8 under a log link) beyond the grid. No random draws.
.landscape_convolve <- function(fitness_model, plan, trait_cols, population_variance,
                                extra_vars, do_clamp, ftype, hull, obs_range) {
    p <- length(trait_cols)
    S <- as.matrix(population_variance)
    # the lattice positions along each trait, counted in steps from the first
    # grid point: only the nodes some grid point's weights reach, and for each
    # grid point (a column of idx) the nodes at offsets -m to m from it
    axes <- lapply(plan$axes, function(a) {
        around <- outer(-a$m:a$m, (seq_len(a$n_grid) - 1) * a$r, "+")
        need <- sort(unique(as.vector(around)))
        c(a, list(need = need, x = a$g1 + need * a$h, idx = matrix(match(around, need), nrow(around)),
                  n_full = (a$n_grid - 1) * a$r + 2 * a$m + 1))
    })
    h <- vapply(axes, function(a) a$h, numeric(1))
    m <- vapply(axes, function(a) a$m, numeric(1))

    # the normal weights on the lattice offsets, correlation included
    offs <- as.matrix(expand.grid(lapply(seq_len(p), function(j) (-m[j]:m[j]) * h[j]), KEEP.OUT.ATTRS = FALSE))
    kern <- exp(-0.5 * rowSums((offs %*% solve(S)) * offs))
    kern <- matrix(kern / sum(kern), 2 * m[1] + 1)

    nodes <- as.matrix(expand.grid(lapply(axes, function(a) a$x), KEEP.OUT.ATTRS = FALSE))
    message(sprintf("Averaging over a lattice of %s nodes, spacing %s",
                    paste(vapply(axes, function(a) length(a$x), integer(1)), collapse = " x "),
                    paste(signif(abs(h), 3), collapse = " x ")))
    # a GAM on the scale of its linear predictor, transformed back
    if (inherits(fitness_model, "gam") && is.function(fitness_model$family$linkinv)) {
        w <- fitness_model$family$linkinv(.predict_fitness(fitness_model, nodes, trait_cols, extra_vars, link = TRUE))
    } else {
        w <- .predict_fitness(fitness_model, nodes, trait_cols, extra_vars)
    }
    clipped <- 0L
    if (do_clamp) {
        held <- .clamp_fitness(w, ftype)
        clipped <- sum(held != w, na.rm = TRUE)
        w <- held
    }

    list(
        mean_fit = .convolve_lattice(w, axes, kern),
        # the share of each population's distribution outside the data
        outside = pmin(pmax(.convolve_lattice(.outside_share(nodes, axes, hull, obs_range), axes, kern), 0), 1),
        clipped = clipped,
        lattice = list(
            n = stats::setNames(vapply(axes, function(a) length(a$x), integer(1)), trait_cols),
            n_full = stats::setNames(vapply(axes, function(a) a$n_full, numeric(1)), trait_cols),
            spacing = stats::setNames(abs(h), trait_cols),
            steps_per_grid = stats::setNames(vapply(axes, function(a) a$r, numeric(1)), trait_cols),
            reach = stats::setNames(as.integer(m), trait_cols),
            reach_sd = plan$reach,
            predictions = nrow(nodes)
        )
    )
}

# the share of each convolve lattice node's cell (the node plus half a step
# either way) lying outside the data: beyond the observed range for one trait,
# worked out exactly, and outside the convex hull of the observed trait pairs
# for two. Taking each node as wholly in or out would be off by up to a few
# hundredths where the hull's edge runs along the lattice, so the cells the
# edge crosses are sampled on a fine grid of points.
.outside_share <- function(nodes, axes, hull, obs_range) {
    hs <- vapply(axes, function(a) abs(a$h), numeric(1))
    if (is.null(hull)) {
        lo <- nodes[, 1] - hs[1] / 2
        hi <- nodes[, 1] + hs[1] / 2
        beyond <- pmax(0, pmin(hi, obs_range[1]) - lo) + pmax(0, hi - pmax(lo, obs_range[2]))
        return(pmin(beyond / hs[1], 1))
    }
    n1 <- length(axes[[1]]$x)
    n2 <- length(axes[[2]]$x)
    inside <- matrix(mgcv::in.out(hull, unname(nodes[, 1:2, drop = FALSE])), n1, n2)
    share <- 1 - inside
    # the cells the edge may cross: a straight edge through a cell parts its
    # node from a neighbouring node along one trait or the other, and a corner
    # of the hull can sit in a cell
    cross <- matrix(FALSE, n1, n2)
    if (n1 > 1L) {
        d <- (inside[-1, , drop = FALSE] != inside[-n1, , drop = FALSE]) & (diff(axes[[1]]$need) == 1)
        cross[-1, ] <- cross[-1, , drop = FALSE] | d
        cross[-n1, ] <- cross[-n1, , drop = FALSE] | d
    }
    if (n2 > 1L) {
        d <- (inside[, -1, drop = FALSE] != inside[, -n2, drop = FALSE]) & rep(diff(axes[[2]]$need) == 1, each = n1)
        cross[, -1] <- cross[, -1, drop = FALSE] | d
        cross[, -n2] <- cross[, -n2, drop = FALSE] | d
    }
    corner <- cbind(match(round((hull[, 1] - axes[[1]]$g1) / axes[[1]]$h), axes[[1]]$need),
                    match(round((hull[, 2] - axes[[2]]$g1) / axes[[2]]$h), axes[[2]]$need))
    cross[corner[stats::complete.cases(corner), , drop = FALSE]] <- TRUE
    e <- which(cross)
    if (length(e)) {
        k <- 16L
        sub <- (seq_len(k) - 0.5) / k - 0.5
        pts <- cbind(rep(nodes[e, 1], each = k * k) + rep(sub * hs[1], times = k),
                     rep(nodes[e, 2], each = k * k) + rep(sub * hs[2], each = k))
        share[e] <- colMeans(matrix(!mgcv::in.out(hull, pts), k * k))
    }
    as.vector(share)
}

# the weighted sum of the lattice values `v` around each grid point, the
# weights `kern` being the normal density on the lattice offsets. One value
# per grid point, the first trait varying fastest as in the grid.
.convolve_lattice <- function(v, axes, kern) {
    a1 <- axes[[1]]
    # for each grid point along the first trait, the lattice nodes its weights reach
    rows <- rep(seq_len(a1$n_grid), each = 2 * a1$m + 1)
    cols <- as.vector(a1$idx)
    if (length(axes) == 1L) {
        return(as.vector(rowsum(rep(kern[, 1], a1$n_grid) * v[cols], rows, reorder = FALSE)))
    }
    a2 <- axes[[2]]
    V <- matrix(v, length(a1$x))
    acc <- matrix(0, a1$n_grid, a2$n_grid)
    # one slice of offsets along the second trait at a time, as a matrix product
    for (b in seq_len(2 * a2$m + 1)) {
        C <- matrix(0, a1$n_grid, length(a1$x))
        C[cbind(rows, cols)] <- rep(kern[, b], a1$n_grid)
        acc <- acc + C %*% V[, a2$idx[b, ], drop = FALSE]
    }
    as.vector(acc)
}
