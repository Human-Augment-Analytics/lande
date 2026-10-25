# Plot Enhanced Correlated Fitness Surface

Plot Enhanced Correlated Fitness Surface

## Usage

``` r
plot_correlated_fitness_enhanced(
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
)
```

## Arguments

- tps:

  Output list from
  [`correlated_fitness_surface()`](https://human-augment-analytics.github.io/lande/dev/reference/correlated_fitness_surface.md).

- trait_cols:

  Character vector of length 2. If `NULL`, inferred from grid.

- original_data:

  Optional data frame of original data.

- fitness_col:

  Optional character string specifying the fitness column for coloring
  points.

- bins:

  Integer specifying the number of contour bins. Default is 12.

- point_alpha:

  Numeric value for point transparency. Default is 0.7.

- show_optimum:

  Logical; mark the highest kept cell, a gold diamond, or an open one
  when it is on the edge of the data. Default is `TRUE`.

- show_groups:

  Logical; when the surface was fitted with a `group`, draw each group's
  mean and the highest point of the surface within that group's hull,
  filled when it is a peak of the surface and open when it is not.
  Default is `TRUE`.

- group_lines:

  Logical; join each group's mean to its peak with a dashed line.
  Default is `TRUE`.

- uncertainty:

  What to draw of the surface's standard error, which a GAM surface
  carries: `"none"` (the default), `"se"` for dashed contour lines of
  the standard error of the fitted fitness over the surface, or `"band"`
  for three panels, the lower bound, the fit and the upper bound, on one
  fill scale.

- ...:

  Additional arguments passed to
  [`ggplot2::labs()`](https://ggplot2.tidyverse.org/reference/labs.html).

## Value

A `ggplot` object with enhanced visualizations.

## Examples

``` r
prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
surf <- correlated_fitness_surface(prep, "survival", c("total_length", "weight"), grid_n = 30)
#> Data type: binary; method: gam; n = 136; k = 29
#> GAM fitting with 136 observations
#>   Trying formula: main
#> Success with formula: main
#> Predictions range: 0.0472 to 0.703
#> Masked 424 of 900 grid points outside the data
plot_correlated_fitness_enhanced(surf, c("total_length", "weight"),
                                 original_data = prep, fitness_col = "survival")
#> Using provided traits: total_length, weight
```
