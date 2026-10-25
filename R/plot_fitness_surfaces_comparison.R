# ======================================================
# plot_fitness_surfaces_comparison.R
# Compare Correlated Fitness Surface vs Adaptive Landscape
#
# IMPORTANT CONCEPT:
# - Correlated Fitness Surface: individual-level (blue)
# - Adaptive Landscape: population-level (red)
# ======================================================

#' Plot Fitness Surfaces Comparison
#'
#' Visualizes the comparison between individual-level correlated fitness 
#' and population-level adaptive landscape.
#'
#' @param comparison_data Output list from \code{compare_fitness_surfaces_data()}.
#' @param bins Integer specifying the number of contour bins. Default is 10.
#' @param title Optional character string for the overarching plot title.
#' @param point_alpha Numeric value for optimum point transparency. Default is 0.7.
#' @param linewidth Numeric value for contour line width. Default is 0.8.
#' @param show_optima Logical indicating whether to show optimum points on the overlay. Default is \code{TRUE}.
#'
#' @return A list containing \code{side_by_side} (a \code{patchwork} object) and \code{overlay} (a \code{ggplot} object).
#' @examples
#' prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
#' surf <- correlated_fitness_surface(prep, "survival", c("total_length", "weight"), grid_n = 15)
#' land <- adaptive_landscape(prep, surf$model, c("total_length", "weight"),
#'                            simulation_n = 100, grid_n = 15)
#' comp <- compare_fitness_surfaces_data(surf, land, c("total_length", "weight"))
#' plot_fitness_surfaces_comparison(comp)
#' @export
plot_fitness_surfaces_comparison <- function(
  comparison_data,
  bins = 10,
  title = NULL,
  point_alpha = 0.7,
  linewidth = 0.8,
  show_optima = TRUE
) {
    # Input validation
    if (!"combined_data" %in% names(comparison_data)) {
        stop("comparison_data must have 'combined_data' element")
    }
    if (!"trait_cols" %in% names(comparison_data)) {
        stop("comparison_data must have 'trait_cols' element")
    }

    combined <- comparison_data$combined_data
    trait_cols <- comparison_data$trait_cols
    optimum_individual <- comparison_data$optimum_individual
    optimum_population <- comparison_data$optimum_population

    if (is.null(title)) {
        title <- "Comparison: Correlated Fitness vs Adaptive Landscape"
    }

    # Split data
    cor_df <- combined[combined$type == "Correlated Fitness (Individual)", ]
    ada_df <- combined[combined$type == "Adaptive Landscape (Population)", ]

    # Left plot: Correlated Fitness (Blue)
    p_cor <- ggplot2::ggplot(cor_df, ggplot2::aes(
        x = .data[[trait_cols[1]]],
        y = .data[[trait_cols[2]]],
        z = .data[["fitness"]]
    )) +
        ggplot2::geom_contour_filled(bins = bins) +
        ggplot2::scale_fill_manual(
            # geom_contour_filled can produce up to bins + 1 bands; supply a
            # couple of extra colours (surplus are ignored) so the fill scale
            # never runs short.
            values = colorRampPalette(c("lightblue", "steelblue", "darkblue", "navy"))(bins + 2),
            name = "Fitness"
        ) +
        ggplot2::labs(
            x = trait_cols[1],
            y = trait_cols[2],
            title = "Correlated Fitness",
            subtitle = "Individual level"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.background = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
            axis.text = ggplot2::element_text(color = "black", size = 10),
            axis.title = ggplot2::element_text(size = 12),
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, color = "gray40"),
            aspect.ratio = 0.8
        )

    # Right plot: Adaptive Landscape (Red)
    p_ada <- ggplot2::ggplot(ada_df, ggplot2::aes(
        x = .data[[trait_cols[1]]],
        y = .data[[trait_cols[2]]],
        z = .data[["fitness"]]
    )) +
        ggplot2::geom_contour_filled(bins = bins) +
        ggplot2::scale_fill_manual(
            # See the note on the correlated-fitness panel: provide a couple of
            # extra colours so the fill scale never runs short.
            values = colorRampPalette(c("lightcoral", "coral", "darkred", "brown"))(bins + 2),
            name = "Mean Fitness"
        ) +
        ggplot2::labs(
            x = trait_cols[1],
            y = trait_cols[2],
            title = "Adaptive Landscape",
            subtitle = "Population level"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.background = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
            axis.text = ggplot2::element_text(color = "black", size = 10),
            axis.title = ggplot2::element_text(size = 12),
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, color = "gray40"),
            aspect.ratio = 0.8
        )

    # Combine side by side with patchwork
    if (requireNamespace("patchwork", quietly = TRUE)) {
        p_side <- p_cor + p_ada +
            patchwork::plot_annotation(
                title = title,
                theme = ggplot2::theme(
                    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14)
                )
            )
    } else {
        warning("Package 'patchwork' not installed. Cannot combine plots.")
        p_side <- NULL
    }

    # Overlay comparison
    p_overlay <- ggplot2::ggplot() +
        # Correlated fitness (blue dashed)
        ggplot2::geom_contour(
            data = cor_df,
            ggplot2::aes(
                x = .data[[trait_cols[1]]],
                y = .data[[trait_cols[2]]],
                z = .data[["fitness"]],
                color = "Correlated Fitness"
            ),
            linetype = "dashed",
            bins = bins,
            linewidth = linewidth
        ) +
        # Adaptive landscape (red solid)
        ggplot2::geom_contour(
            data = ada_df,
            ggplot2::aes(
                x = .data[[trait_cols[1]]],
                y = .data[[trait_cols[2]]],
                z = .data[["fitness"]],
                color = "Adaptive Landscape"
            ),
            linetype = "solid",
            bins = bins,
            linewidth = linewidth
        )

    # Optimum points, only when requested and available
    optimum_layer <- function(opt, label) {
        if (!show_optima || is.null(opt) || !all(trait_cols %in% names(opt))) {
            return(NULL)
        }
        ggplot2::geom_point(
            data = opt,
            ggplot2::aes(
                x = .data[[trait_cols[1]]],
                y = .data[[trait_cols[2]]],
                color = label
            ),
            size = 4,
            shape = 18,
            alpha = point_alpha
        )
    }
    p_overlay <- p_overlay +
        optimum_layer(optimum_individual, "Individual Optimum") +
        optimum_layer(optimum_population, "Population Optimum") +
        # Color scale
        ggplot2::scale_color_manual(
            name = "",
            values = c(
                "Correlated Fitness" = "steelblue",
                "Adaptive Landscape" = "darkred",
                "Individual Optimum" = "darkblue",
                "Population Optimum" = "brown"
            ),
            labels = c(
                "Correlated Fitness" = "Correlated Fitness (Individual)",
                "Adaptive Landscape" = "Adaptive Landscape (Population)",
                "Individual Optimum" = "Individual Optimum",
                "Population Optimum" = "Population Optimum"
            )
        ) +
        # Labels
        ggplot2::labs(
            x = trait_cols[1],
            y = trait_cols[2],
            title = "Overlay Comparison",
            subtitle = "Blue dashed: Individual fitness | Red solid: Population mean fitness\nBlue star: Individual optimum | Red star: Population optimum"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.background = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
            axis.text = ggplot2::element_text(color = "black", size = 10),
            axis.title = ggplot2::element_text(size = 12),
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, color = "gray40"),
            legend.position = "right",
            legend.box = "vertical",
            aspect.ratio = 1.0
        )

    return(list(
        side_by_side = p_side,
        overlay = p_overlay
    ))
}
