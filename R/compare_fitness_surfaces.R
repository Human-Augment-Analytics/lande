# ======================================================
# compare_fitness_surfaces_data.R
# Compare Correlated Fitness Surface vs Adaptive Landscape
#
# IMPORTANT CONCEPT:
# - Correlated Fitness Surface: individual-level fitness (w ~ z)
# - Adaptive Landscape: population-level mean fitness (Wbar ~ zbar)
#
# KEY PRINCIPLE:
# - Both surfaces should be based on standardized traits
# - Use prepare_selection_data() before creating both surfaces
# - DO NOT standardize again within this function
# ======================================================

#' Compare Correlated Fitness Surface vs Adaptive Landscape
#'
#' Compares individual-level fitness (correlated surface) to population-level mean fitness (adaptive landscape).
#'
#' @param correlated_surface Output list from \code{correlated_fitness_surface()}.
#' @param adaptive_landscape Output list from adaptive landscape modeling.
#' @param trait_cols A character vector of length 2 specifying the trait column names.
#' @param calculate_correlation Logical indicating whether to calculate the correlation between the two surfaces.
#'
#' @return A list containing combined grid data, summary stats, optima points, and correlation metrics.
#' @export
#'
#' @examples
#' prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
#' surf <- correlated_fitness_surface(prep, "survival", c("total_length", "weight"), grid_n = 15)
#' land <- adaptive_landscape(prep, surf$model, c("total_length", "weight"),
#'                            simulation_n = 100, grid_n = 15)
#' comp <- compare_fitness_surfaces_data(surf, land, c("total_length", "weight"))
#' plot_fitness_surfaces_comparison(comp)
compare_fitness_surfaces_data <- function(
  correlated_surface,
  adaptive_landscape,
  trait_cols,
  calculate_correlation = TRUE
) {
    # Input validation
    if (!"grid" %in% names(correlated_surface)) {
        stop("correlated_surface must have a 'grid' element")
    }
    if (!"grid" %in% names(adaptive_landscape)) {
        stop("adaptive_landscape must have a 'grid' element")
    }

    cor_df <- correlated_surface$grid

    # Check if required columns exist (trait columns should already be named)
    if (!all(trait_cols %in% names(cor_df))) {
        stop(
            "Trait columns not found in correlated_surface$grid. ",
            "Expected: ", paste(trait_cols, collapse = ", "), "\n",
            "Found: ", paste(names(cor_df), collapse = ", ")
        )
    }
    if (!".fit" %in% names(cor_df)) {
        stop("No '.fit' column found in correlated_surface$grid")
    }

    # Add fitness and type columns (preserve original trait column names)
    cor_df$fitness <- cor_df$.fit
    cor_df$type <- "Correlated Fitness (Individual)"

    ada_df <- adaptive_landscape$grid

    if (!all(trait_cols %in% names(ada_df))) {
        stop(
            "Trait columns not found in adaptive_landscape$grid. ",
            "Expected: ", paste(trait_cols, collapse = ", "), "\n",
            "Found: ", paste(names(ada_df), collapse = ", ")
        )
    }
    if (!".mean_fit" %in% names(ada_df)) {
        stop("No '.mean_fit' column found in adaptive_landscape$grid")
    }

    ada_df$fitness <- ada_df$.mean_fit
    ada_df$type <- "Adaptive Landscape (Population)"

    n_cor <- nrow(cor_df)
    n_ada <- nrow(ada_df)
    if (n_cor != n_ada) {
        warning(
            "Grid sizes differ: correlated surface has ", n_cor,
            " points, adaptive landscape has ", n_ada, " points.\n",
            "Comparison may be affected. Consider using same grid_n."
        )
    }

    if (requireNamespace("dplyr", quietly = TRUE)) {
        combined <- dplyr::bind_rows(
            cor_df[, c(trait_cols, "fitness", "type")],
            ada_df[, c(trait_cols, "fitness", "type")]
        )
    } else {
        # Fallback to rbind (risk of factor issues)
        warning("dplyr not available. Using rbind - may cause issues with factors.")
        combined <- rbind(
            cor_df[, c(trait_cols, "fitness", "type")],
            ada_df[, c(trait_cols, "fitness", "type")]
        )
    }

    summary_stats <- data.frame(
        Surface = c("Correlated Fitness", "Adaptive Landscape"),
        Fitness_Range_Min = c(
            min(cor_df$fitness, na.rm = TRUE),
            min(ada_df$fitness, na.rm = TRUE)
        ),
        Fitness_Range_Max = c(
            max(cor_df$fitness, na.rm = TRUE),
            max(ada_df$fitness, na.rm = TRUE)
        ),
        Fitness_Mean = c(
            mean(cor_df$fitness, na.rm = TRUE),
            mean(ada_df$fitness, na.rm = TRUE)
        ),
        Fitness_SD = c(
            sd(cor_df$fitness, na.rm = TRUE),
            sd(ada_df$fitness, na.rm = TRUE)
        ),
        N_Points = c(n_cor, n_ada)
    )

    # Find optimum points
    optimum_individual <- cor_df[which.max(cor_df$fitness), ]
    optimum_population <- ada_df[which.max(ada_df$fitness), ]

    dist_opt <- sqrt(
        sum((optimum_individual[, trait_cols] - optimum_population[, trait_cols])^2)
    )

    correlation <- NULL
    if (calculate_correlation && n_cor == n_ada) {
        # Need to align points by trait values
        if (requireNamespace("akima", quietly = TRUE)) {
            # Create common grid
            x_range <- range(c(cor_df[[trait_cols[1]]], ada_df[[trait_cols[1]]]))
            y_range <- range(c(cor_df[[trait_cols[2]]], ada_df[[trait_cols[2]]]))

            x_grid <- seq(x_range[1], x_range[2], length.out = 50)
            y_grid <- seq(y_range[1], y_range[2], length.out = 50)

            # Interpolate both surfaces (masked surface cells are NA and are
            # left out; interp works from the remaining irregular points)
            cor_ok <- cor_df[!is.na(cor_df$fitness), , drop = FALSE]
            cor_interp <- akima::interp(
                cor_ok[[trait_cols[1]]], cor_ok[[trait_cols[2]]], cor_ok$fitness,
                xo = x_grid, yo = y_grid, linear = TRUE
            )
            ada_interp <- akima::interp(
                ada_df[[trait_cols[1]]], ada_df[[trait_cols[2]]], ada_df$fitness,
                xo = x_grid, yo = y_grid, linear = TRUE
            )

            correlation <- cor(as.vector(cor_interp$z), as.vector(ada_interp$z),
                use = "complete.obs"
            )
        } else {
            # Fallback: use raw grid points (may not be aligned)
            warning("Package 'akima' not installed. Correlation may be approximate.")
            correlation <- cor(cor_df$fitness, ada_df$fitness, use = "complete.obs")
        }
    }

    .msg_table("Individual optimum:", optimum_individual[, c(trait_cols, "fitness")])
    .msg_table("Population optimum:", optimum_population[, c(trait_cols, "fitness")])
    message("Distance between optima: ", round(dist_opt, 4))
    .msg_table("Summary statistics:", summary_stats)
    if (!is.null(correlation)) {
        message("Correlation between surfaces: ", round(correlation, 4))
    }


    list(
        combined_data = combined,
        summary_stats = summary_stats,
        optimum_individual = optimum_individual,
        optimum_population = optimum_population,
        distance_between_optima = dist_opt,
        correlation = correlation,
        correlated_surface = correlated_surface,
        adaptive_landscape = adaptive_landscape,
        trait_cols = trait_cols,
        grid_aligned = (n_cor == n_ada)
    )
}
