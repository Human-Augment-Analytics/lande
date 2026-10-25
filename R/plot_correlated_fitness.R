# ======================================================
# plot_correlated_fitness.R
# Visualize correlated fitness surface (individual-level)
#
# IMPORTANT CONCEPT:
# - Correlated Fitness Surface: individual fitness (w ~ z1, z2)
# - This is DIFFERENT from adaptive landscape (population-level)
#
# KEY PRINCIPLE:
# - Traits should already be standardized
# - DO NOT standardize again within this function
# ======================================================

#' @noRd
# internal utility: what to draw for a surface. When the grid carries the
# unmasked predictions the full surface is drawn and the blanked region is
# covered: the outside of the hull by a polygon, cells too far from any
# individual by a white contour band on the distance field. Both give clean
# edges. Otherwise masked (NA) cells are simply dropped.
.surface_layers <- function(tps, df, trait1, trait2, fit_col) {
  if (!".fit_all" %in% names(df)) {
    return(list(df = df[!is.na(df[[fit_col]]), , drop = FALSE], fit_col = fit_col, cover = NULL))
  }
  cover <- list()
  hull <- tps$hull
  if (!is.null(hull)) {
    frame <- .outside_frame(hull, trait1, trait2, range(df[[trait1]]), range(df[[trait2]]))
    cover <- c(cover, list(ggplot2::geom_polygon(
      data = frame,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill = "white", colour = "white", linewidth = 0.5, inherit.aes = FALSE
    )))
  }
  too_far <- tps$too_far
  if (!is.null(too_far) && ".dist" %in% names(df) && any(df$.dist > too_far)) {
    # two bands of the distance field, near and far; only the far one is drawn
    cover <- c(cover, list(
      ggplot2::layer(
        stat = .StatDistanceBand, geom = "polygon", position = "identity",
        data = df,
        mapping = ggplot2::aes(
          x = .data[[trait1]], y = .data[[trait2]], z = .data$.dist,
          alpha = ggplot2::after_stat(.data$level)
        ),
        params = list(breaks = c(-1, too_far, max(df$.dist) + 1), fill = "white", na.rm = FALSE),
        inherit.aes = FALSE
      ),
      ggplot2::scale_alpha_manual(values = c(0, 1), guide = "none")
    ))
  }
  list(df = df, fit_col = ".fit_all", cover = cover)
}

#' @noRd
# internal utility: the region outside the hull as one ring: round the grid
# rectangle, along a bridge to the nearest hull vertex, once round the hull
# the opposite way, and back. Renders as a frame with a hull-shaped window.
.outside_frame <- function(hull, trait1, trait2, xr, yr) {
  h <- cbind(hull[[trait1]], hull[[trait2]])
  h <- h[-nrow(h), , drop = FALSE]
  corner <- c(xr[1], yr[1])
  i <- which.min((h[, 1] - corner[1])^2 + (h[, 2] - corner[2])^2)
  ring <- h[c(i:nrow(h), seq_len(i)), , drop = FALSE]
  frame <- rbind(
    c(xr[1], yr[1]), c(xr[2], yr[1]), c(xr[2], yr[2]), c(xr[1], yr[2]), c(xr[1], yr[1]),
    ring, c(xr[1], yr[1])
  )
  data.frame(x = frame[, 1], y = frame[, 2])
}

#' @noRd
# the filled-contour stat without its fill mapping, so the white band drawn
# over cells far from the data adds no keys to the fitness legend
.StatDistanceBand <- ggplot2::ggproto(
  "StatDistanceBand", ggplot2::StatContourFilled,
  default_aes = ggplot2::aes(order = ggplot2::after_stat(level))
)

#' @noRd
# internal utility: what the uncertainty option adds to a surface plot. "se"
# overlays dashed contour lines of the standard error of the fitted fitness;
# "band" redraws the surface as three panels, lower, fit and upper, on one
# fill scale. A surface without standard errors is drawn as the fit alone.
.uncertainty <- function(tps, df, z_col, trait1, trait2, uncertainty) {
  out <- list(df = df, z_col = z_col, se_layer = NULL, extra = NULL)
  if (uncertainty == "none") return(out)
  if (!".se" %in% names(tps$grid) || all(is.na(tps$grid$.se))) {
    warning("This surface has no standard errors (they come with method = \"gam\"); drawing the fit alone")
    return(out)
  }
  full <- z_col == ".fit_all"
  if (uncertainty == "se") {
    se_col <- if (full) ".se_all" else ".se"
    out$se_layer <- ggplot2::geom_contour(
      data = df, ggplot2::aes(x = .data[[trait1]], y = .data[[trait2]], z = .data[[se_col]]),
      colour = "white", linetype = "dashed", linewidth = 0.4, bins = 6, inherit.aes = FALSE
    )
    out$extra <- ggplot2::labs(caption = "Dashed white: standard error of the fitted fitness")
    return(out)
  }
  cols <- c(if (full) ".fit_lo_all" else ".fit_lo", z_col, if (full) ".fit_hi_all" else ".fit_hi")
  panels <- c("Lower", "Fit", "Upper")
  out$df <- do.call(rbind, lapply(1:3, function(i) {
    d <- df
    d$.z <- d[[cols[i]]]
    d$.panel <- factor(panels[i], levels = panels)
    d
  }))
  out$z_col <- ".z"
  out$extra <- ggplot2::facet_wrap(~ .panel)
  out
}

#' @noRd
# internal utility: gold diamond for the optimum, open if it is on the edge
.optimum_layer <- function(tps, df, trait1, trait2, fit_col) {
  col <- if (".fit" %in% names(df)) ".fit" else fit_col
  d <- df[!is.na(df[[col]]), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  opt <- d[which.max(d[[col]]), , drop = FALSE]
  pk <- tps$peaks
  interior <- if (is.null(pk) || !nrow(pk)) TRUE else {
    same <- abs(pk[[trait1]] - opt[[trait1]]) < 1e-8 & abs(pk[[trait2]] - opt[[trait2]]) < 1e-8
    any(same & pk$interior)
  }
  if (interior) {
    ggplot2::geom_point(data = opt, ggplot2::aes(x = .data[[trait1]], y = .data[[trait2]], shape = "Optimum"),
                        colour = "gold", size = 4, inherit.aes = FALSE)
  } else {
    ggplot2::geom_point(data = opt, ggplot2::aes(x = .data[[trait1]], y = .data[[trait2]], shape = "Edge maximum"),
                        colour = "black", size = 3.5, stroke = 0.9, inherit.aes = FALSE)
  }
}

#' @noRd
# internal utility: group means (open circles, labelled) and the highest
# point of the surface within each group's hull, a filled triangle when that
# point is a peak of the surface and an open one when the surface keeps
# rising past the group's range or the edge of the data, joined to the mean
# by a dashed line. Only present when the surface was fitted with a group.
.group_layers <- function(tps, trait1, trait2, lines = TRUE) {
  g <- tps$groups
  if (is.null(g) || !nrow(g)) return(NULL)
  m1 <- paste0("mean_", trait1)
  m2 <- paste0("mean_", trait2)
  p1 <- paste0("peak_", trait1)
  p2 <- paste0("peak_", trait2)
  if (!all(c(m1, m2, p1, p2) %in% names(g))) return(NULL)
  with_peak <- g[!is.na(g[[p1]]), , drop = FALSE]
  # surfaces fitted before the flag existed have no peak_interior column
  is_peak <- if ("peak_interior" %in% names(with_peak)) with_peak$peak_interior %in% TRUE else rep(TRUE, nrow(with_peak))
  peaks <- with_peak[is_peak, , drop = FALSE]
  edges <- with_peak[!is_peak, , drop = FALSE]
  list(
    if (lines) ggplot2::geom_segment(
      data = with_peak,
      ggplot2::aes(x = .data[[m1]], y = .data[[m2]], xend = .data[[p1]], yend = .data[[p2]]),
      colour = "black", linewidth = 0.4, linetype = "dashed", inherit.aes = FALSE
    ),
    ggplot2::geom_point(
      data = g, ggplot2::aes(x = .data[[m1]], y = .data[[m2]], shape = "Group mean"),
      fill = "white", colour = "black", size = 3, inherit.aes = FALSE
    ),
    if (nrow(peaks)) ggplot2::geom_point(
      data = peaks, ggplot2::aes(x = .data[[p1]], y = .data[[p2]], shape = "Group peak"),
      fill = "black", colour = "white", size = 2.8, inherit.aes = FALSE
    ),
    if (nrow(edges)) ggplot2::geom_point(
      data = edges, ggplot2::aes(x = .data[[p1]], y = .data[[p2]], shape = "Group high point"),
      colour = "black", size = 2.6, stroke = 0.9, inherit.aes = FALSE
    ),
    ggplot2::geom_text(
      data = .shared_labels(g, c(m1, m2), "group"), ggplot2::aes(x = .data[[m1]], y = .data[[m2]], label = .data$group),
      vjust = -1, size = 3, colour = "black", inherit.aes = FALSE
    )
  )
}

#' @noRd
# internal utility: one label per position, so groups that share a mean get
# "a, b" instead of labels printed over each other. With each group
# standardised within itself every mean is zero.
.shared_labels <- function(df, cols, label) {
  key <- do.call(paste, lapply(cols, function(k) round(df[[k]], 6)))
  first <- !duplicated(key)
  out <- df[first, , drop = FALSE]
  out[[label]] <- vapply(
    key[first],
    function(k) paste(as.character(df[[label]][key == k]), collapse = ", "),
    character(1), USE.NAMES = FALSE
  )
  rownames(out) <- NULL
  out
}

#' @noRd
# one key for the marks drawn on a surface: the optimum, or the highest cell
# when it sits at the edge, and, with a group, the group means and peaks.
# Only the marks present appear in it.
.mark_key <- function() {
  ggplot2::scale_shape_manual(
    name = NULL,
    values = c("Optimum" = 18, "Edge maximum" = 5, "Group mean" = 21, "Group peak" = 24, "Group high point" = 2),
    guide = ggplot2::guide_legend(order = 1)
  )
}

#' Plot Correlated Fitness Surface
#'
#' @param tps Output list from \code{correlated_fitness_surface()}.
#' @param trait_cols Character vector of length 2 specifying the trait column names.
#' @param bins Integer specifying the number of contour bins. Default is 12.
#' @param point_alpha Numeric value for point transparency. Default is 0.7.
#' @param show_points Logical indicating whether to show original data points. Default is \code{FALSE}.
#' @param show_optimum Logical; mark the highest kept cell, a gold diamond, or
#'   an open one when it is on the edge of the data. Default is \code{TRUE}.
#' @param show_groups Logical; when the surface was fitted with a \code{group},
#'   draw each group's mean (open circle, labelled) and the highest point of
#'   the surface within that group's hull: a filled triangle when it is a
#'   peak of the surface, an open one when the surface keeps rising past the
#'   group's range or the edge of the data. Default is \code{TRUE}.
#' @param group_lines Logical; join each group's mean to its peak with a dashed
#'   line. Default is \code{TRUE}.
#' @param uncertainty What to draw of the surface's standard error, which a GAM
#'   surface carries: \code{"none"} (the default), \code{"se"} for dashed
#'   contour lines of the standard error of the fitted fitness over the
#'   surface, or \code{"band"} for three panels, the lower bound, the fit and
#'   the upper bound, on one fill scale.
#' @param ... Additional arguments passed to \code{ggplot2::labs()}.
#'
#' @return A \code{ggplot} object representing the correlated fitness surface.
#' @examples
#' prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
#' surf <- correlated_fitness_surface(prep, "survival", c("total_length", "weight"), grid_n = 30)
#' plot_correlated_fitness(surf, c("total_length", "weight"), show_points = TRUE)
#' @export
plot_correlated_fitness <- function(
  tps,
  trait_cols,
  bins = 12,
  point_alpha = 0.7,
  show_points = FALSE,
  show_optimum = TRUE,
  show_groups = TRUE,
  group_lines = TRUE,
  uncertainty = c("none", "se", "band"),
  ...
) {
  # Input validation
  stopifnot(is.list(tps), "grid" %in% names(tps))
  uncertainty <- match.arg(uncertainty)

  df <- tps$grid

  if (length(trait_cols) != 2) {
    stop("trait_cols must be length 2.")
  }
  if (!all(trait_cols %in% names(df))) {
    stop("Trait columns not found in tps$grid.")
  }

  # Find fitness column
  if (".fit" %in% names(df)) {
    fitness_col <- ".fit"
  } else if ("fitness" %in% names(df)) {
    fitness_col <- "fitness"
  } else if ("pred" %in% names(df)) {
    fitness_col <- "pred"
  } else {
    stop("No fitness column found in grid.")
  }

  # Masked surfaces: draw the full surface, then cover the outside of the data hull
  layers <- .surface_layers(tps, df, trait_cols[1], trait_cols[2], fitness_col)
  unc <- .uncertainty(tps, layers$df, layers$fit_col, trait_cols[1], trait_cols[2], uncertainty)
  draw_df <- unc$df
  z_col <- unc$z_col

  p <- ggplot2::ggplot(draw_df) +
    # Filled contours
    ggplot2::geom_contour_filled(
      ggplot2::aes(
        x = .data[[trait_cols[1]]],
        y = .data[[trait_cols[2]]],
        z = .data[[z_col]]
      ),
      bins = bins,
      inherit.aes = FALSE
    ) +
    # Contour lines
    ggplot2::geom_contour(
      ggplot2::aes(
        x = .data[[trait_cols[1]]],
        y = .data[[trait_cols[2]]],
        z = .data[[z_col]]
      ),
      color = "black",
      alpha = 0.3,
      linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    unc$se_layer +
    layers$cover +
    # Labels
    ggplot2::labs(
      x = trait_cols[1],
      y = trait_cols[2],
      fill = "Fitness",
      title = "Correlated Fitness Surface",
      subtitle = "Individual-level fitness",
      ...
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
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, color = "gray40")
    )

  # Add data points if requested
  if (show_points && "original_data" %in% names(tps)) {
    orig_data <- tps$original_data
    if (!is.null(orig_data) && all(trait_cols %in% names(orig_data))) {
      p <- p + ggplot2::geom_point(
        data = orig_data,
        ggplot2::aes(
          x = .data[[trait_cols[1]]],
          y = .data[[trait_cols[2]]]
        ),
        size = 1.5,
        alpha = point_alpha,
        color = "gray40"
      )
    }
  }

  if (!is.null(unc$extra)) p <- p + unc$extra

  if (show_optimum) {
    p <- p + .optimum_layer(tps, df, trait_cols[1], trait_cols[2], fitness_col)
  }

  if (show_groups) {
    p <- p + .group_layers(tps, trait_cols[1], trait_cols[2], lines = group_lines)
  }
  p <- p + .mark_key()

  return(p)
}


# ======================================================
# plot_correlated_fitness_enhanced
# Enhanced version with original data points colored by fitness
# ======================================================

#' Plot Enhanced Correlated Fitness Surface
#'
#' @param tps Output list from \code{correlated_fitness_surface()}.
#' @param trait_cols Character vector of length 2. If \code{NULL}, inferred from grid.
#' @param original_data Optional data frame of original data.
#' @param fitness_col Optional character string specifying the fitness column for coloring points.
#' @param bins Integer specifying the number of contour bins. Default is 12.
#' @param point_alpha Numeric value for point transparency. Default is 0.7.
#' @param show_groups Logical; when the surface was fitted with a \code{group},
#'   draw each group's mean and the highest point of the surface within that
#'   group's hull, filled when it is a peak of the surface and open when it
#'   is not. Default is \code{TRUE}.
#' @param ... Additional arguments passed to \code{ggplot2::labs()}.
#' @inheritParams plot_correlated_fitness
#'
#' @return A \code{ggplot} object with enhanced visualizations.
#' @examples
#' prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
#' surf <- correlated_fitness_surface(prep, "survival", c("total_length", "weight"), grid_n = 30)
#' plot_correlated_fitness_enhanced(surf, c("total_length", "weight"),
#'                                  original_data = prep, fitness_col = "survival")
#' @export
plot_correlated_fitness_enhanced <- function(
  tps,
  trait_cols = NULL,
  original_data = NULL,
  fitness_col = NULL,
  bins = 12,
  point_alpha = 0.7,
  show_optimum = TRUE,
  show_groups = TRUE,
  group_lines = TRUE,
  uncertainty = c("none", "se", "band"),
  ...
) {
  # Input validation
  uncertainty <- match.arg(uncertainty)
  df <- tps$grid

  # Determine trait columns
  if (!is.null(trait_cols) && length(trait_cols) == 2) {
    trait1 <- trait_cols[1]
    trait2 <- trait_cols[2]
    message("Using provided traits: ", trait1, ", ", trait2)
  } else if (!is.null(tps$trait_cols) && length(tps$trait_cols) == 2) {
    trait1 <- tps$trait_cols[1]
    trait2 <- tps$trait_cols[2]
    message("Using traits from surface object: ", trait1, ", ", trait2)
  } else {
    # Infer from grid columns, excluding any grouping column carried along
    possible_traits <- setdiff(
      names(df),
      c(".fit", ".fit_all", ".inside", ".dist", ".se", ".se_all", ".fit_lo", ".fit_lo_all", ".fit_hi", ".fit_hi_all",
        "fitness", "pred", "fit", "lwr", "upr", "type", "surface_type", tps$group_used)
    )
    if (length(possible_traits) >= 2) {
      trait1 <- possible_traits[1]
      trait2 <- possible_traits[2]
      message("Inferred traits: ", trait1, ", ", trait2)
    } else {
      stop("Cannot determine trait columns. Please provide trait_cols.")
    }
  }

  # Determine fitness column in grid
  if (".fit" %in% names(df)) {
    fit_col <- ".fit"
  } else if ("fitness" %in% names(df)) {
    fit_col <- "fitness"
  } else {
    stop("No fitness column found in grid")
  }

  # Masked surfaces: draw the full surface, then cover the outside of the data hull
  layers <- .surface_layers(tps, df, trait1, trait2, fit_col)
  unc <- .uncertainty(tps, layers$df, layers$fit_col, trait1, trait2, uncertainty)
  draw_df <- unc$df
  z_col <- unc$z_col

  p <- ggplot2::ggplot() +
    ggplot2::geom_contour_filled(
      data = draw_df,
      ggplot2::aes(
        x = .data[[trait1]],
        y = .data[[trait2]],
        z = .data[[z_col]]
      ),
      bins = bins
    ) +
    ggplot2::geom_contour(
      data = draw_df,
      ggplot2::aes(
        x = .data[[trait1]],
        y = .data[[trait2]],
        z = .data[[z_col]]
      ),
      color = "black",
      alpha = 0.3,
      linewidth = 0.3
    ) +
    unc$se_layer +
    layers$cover +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.text = ggplot2::element_text(color = "black", size = 10),
      axis.title = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, color = "gray40")
    )

  if (!is.null(original_data) && !is.null(fitness_col)) {
    if (all(c(trait1, trait2, fitness_col) %in% names(original_data))) {
      unique_vals <- unique(original_data[[fitness_col]])
      is_binary <- length(unique_vals) == 2 && all(unique_vals %in% c(0, 1))

      if (is_binary) {
        p <- p +
          ggplot2::geom_point(
            data = original_data,
            ggplot2::aes(
              x = .data[[trait1]],
              y = .data[[trait2]],
              color = as.factor(.data[[fitness_col]])
            ),
            size = 2,
            alpha = point_alpha
          ) +
          ggplot2::scale_color_manual(
            values = c("0" = "#D55E00", "1" = "#0072B2"),
            labels = c("0" = "Perished", "1" = "Survived"),
            name = "Outcome"
          )
      } else {
        p <- p +
          ggplot2::geom_point(
            data = original_data,
            ggplot2::aes(
              x = .data[[trait1]],
              y = .data[[trait2]],
              color = .data[[fitness_col]]
            ),
            size = 2,
            alpha = point_alpha
          ) +
          ggplot2::scale_color_viridis_c(name = "Fitness")
      }
    } else {
      warning("Required columns not found in original_data")
    }
  }

  if (!is.null(unc$extra)) p <- p + unc$extra

  if (show_optimum) {
    p <- p + .optimum_layer(tps, df, trait1, trait2, fit_col)
  }

  if (show_groups) {
    p <- p + .group_layers(tps, trait1, trait2, lines = group_lines)
  }
  p <- p + .mark_key()

  p <- p + ggplot2::labs(
    x = trait1,
    y = trait2,
    title = "Correlated Fitness Surface",
    subtitle = "Individual-level fitness",
    ...
  )

  return(p)
}


# ======================================================
# plot_correlated_fitness_3d
# 3D perspective view of the correlated fitness surface
# ======================================================

#' Plot Correlated Fitness Surface (3D Perspective)
#'
#' Draws the fitted fitness of a two-trait surface in perspective. Cells
#' blanked by \code{mask} or \code{too_far} are left out, so only the part of
#' the surface with data under it is drawn. Each grid cell is outlined, so a
#' surface fitted with \code{grid_n = 80} shows as 79 x 79 small squares.
#'
#' @param tps Output list from \code{correlated_fitness_surface()}.
#' @param trait_cols Character vector of length 2 specifying the trait column names.
#' @param theta Numeric azimuthal viewing angle. Default is -30.
#' @param phi Numeric colatitude viewing angle. Default is 30.
#' @param color_palette Colours for the surface: a vector of colours, or
#'   \code{"jet"} for the rainbow ramp (dark blue through cyan, yellow and red
#'   to dark red) used in many published fitness surfaces. Defaults to
#'   viridis. Jet is harder to read for colour-blind viewers and its bright
#'   bands can look like steps that are not in the surface, so viridis is the
#'   safer choice for most figures.
#' @param show_peaks Logical; mark the interior peaks of the surface (the
#'   rows of \code{tps$peaks} with \code{interior} \code{TRUE}) with gold
#'   diamonds. Maxima at the edge of the data are not marked. Default is
#'   \code{TRUE}.
#' @param border Colour of the lines outlining each grid cell. Default is
#'   \code{"grey30"}, a dark grey that keeps the squares clear without
#'   darkening the colours as black lines do; \code{NA} draws the surface
#'   without them. With more than 100 grid values along a side the default
#'   lines are faded, more as the grid gets finer, so they do not hide the
#'   colours; give \code{border} to keep them at full strength.
#' @param lwd Width of the cell outlines, which also sets the width of the box
#'   lines. Default is 0.5.
#' @param shade Strength of the shading that shows the relief, from 0 (flat
#'   colours, no shading) towards 1 (darker slopes), passed to
#'   \code{graphics::persp()}. Default is 0.15, light enough that the colours
#'   stay close to the colour bar while the hills still stand out. Use
#'   \code{NA} to switch the shading off.
#' @param ... Additional arguments passed to \code{fields::drape.plot()}.
#'
#' @return The viewing transformation matrix from \code{fields::drape.plot()},
#'   invisibly, which \code{grDevices::trans3d()} uses to add points.
#' @examples
#' prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
#' surf <- correlated_fitness_surface(prep, "survival", c("total_length", "weight"), grid_n = 30)
#' plot_correlated_fitness_3d(surf, c("total_length", "weight"))
#' @export
plot_correlated_fitness_3d <- function(
  tps,
  trait_cols,
  theta = -30,
  phi = 30,
  color_palette = NULL,
  show_peaks = TRUE,
  border = "grey30",
  lwd = 0.5,
  shade = 0.15,
  ...
) {
  # Input validation
  stopifnot(is.list(tps), "grid" %in% names(tps))

  if (length(trait_cols) != 2L) {
    stop("The 3D view needs two traits; for one trait use univariate_spline() ",
         "and plot_univariate_fitness()")
  }

  if (!requireNamespace("fields", quietly = TRUE)) {
    stop("Package 'fields' is required for 3D surface plots. ",
         "Install it with install.packages('fields').")
  }

  df <- tps$grid

  if (!all(trait_cols %in% names(df))) {
    stop("Trait columns not found in tps$grid.")
  }
  if (!".fit" %in% names(df)) {
    stop("No fitness column found in grid.")
  }

  x <- df[[trait_cols[1]]]
  y <- df[[trait_cols[2]]]

  # correlated_fitness_surface() predicts on a regular grid, so the fit goes
  # straight into a matrix. Masked cells stay NA and persp() leaves them out.
  xu <- sort(unique(x))
  yu <- sort(unique(y))
  if (nrow(df) != length(xu) * length(yu)) {
    stop("tps$grid is not a complete regular grid; cannot draw a 3D surface")
  }
  zmat <- matrix(NA_real_, nrow = length(xu), ncol = length(yu))
  zmat[cbind(match(x, xu), match(y, yu))] <- df$.fit
  if (all(is.na(zmat))) {
    stop("Every cell of the surface is masked; nothing to draw")
  }

  if (is.null(color_palette)) {
    if (requireNamespace("viridis", quietly = TRUE)) {
      color_palette <- viridis::viridis(100)
    } else {
      color_palette <- grDevices::heat.colors(100)
    }
  } else if (identical(color_palette, "jet")) {
    color_palette <- jet_colours(100)
  }

  if (missing(border)) border <- grid_line_colour(border, max(dim(zmat)))

  pm <- fields::drape.plot(
    x = xu,
    y = yu,
    z = zmat,
    theta = theta,
    phi = phi,
    xlab = trait_cols[1],
    ylab = trait_cols[2],
    zlab = "Fitness",
    main = "Correlated Fitness Surface (3D View)",
    col = color_palette,
    border = border,
    lwd = lwd,
    shade = shade,
    ...
  )

  # interior peaks, the same gold diamond as the optimum in the 2D plots
  pk <- tps$peaks
  if (show_peaks && !is.null(pk) && nrow(pk) && all(c(trait_cols, "fit", "interior") %in% names(pk))) {
    pk <- pk[pk$interior %in% TRUE, , drop = FALSE]
    if (nrow(pk)) {
      at <- grDevices::trans3d(pk[[trait_cols[1]]], pk[[trait_cols[2]]], pk$fit, pm)
      graphics::points(at, pch = 23, bg = "gold", col = "black", cex = 1.4)
    }
  }

  invisible(pm)
}


# Past 100 grid values a side the cell outlines would hide the colours, so the
# default line colour fades in step with the grid.
grid_line_colour <- function(border, n_side) {
  if (is.na(border) || n_side <= 100) return(border)
  grDevices::adjustcolor(border, alpha.f = (100 / n_side)^2)
}


# The classic jet rainbow: dark blue through cyan, yellow and red to dark red.
# fields::tim.colors() has since moved to a turbo-like ramp, so jet is built
# here from its nine anchor colours.
jet_colours <- function(n) {
  grDevices::colorRampPalette(c("#00007F", "#0000FF", "#007FFF", "#00FFFF", "#7FFF7F",
                                "#FFFF00", "#FF7F00", "#FF0000", "#7F0000"))(n)
}
