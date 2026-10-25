# Compare Correlated Fitness Surface vs Adaptive Landscape

Compares individual-level fitness (correlated surface) to
population-level mean fitness (adaptive landscape).

## Usage

``` r
compare_fitness_surfaces_data(
  correlated_surface,
  adaptive_landscape,
  trait_cols,
  calculate_correlation = TRUE
)
```

## Arguments

- correlated_surface:

  Output list from
  [`correlated_fitness_surface()`](https://human-augment-analytics.github.io/lande/dev/reference/correlated_fitness_surface.md).

- adaptive_landscape:

  Output list from adaptive landscape modeling.

- trait_cols:

  A character vector of length 2 specifying the trait column names.

- calculate_correlation:

  Logical indicating whether to calculate the correlation between the
  two surfaces.

## Value

A list containing combined grid data, summary stats, optima points, and
correlation metrics.

## Examples

``` r
prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
surf <- correlated_fitness_surface(prep, "survival", c("total_length", "weight"), grid_n = 15)
#> Data type: binary; method: gam; n = 136; k = 29
#> GAM fitting with 136 observations
#>   Trying formula: main
#> Success with formula: main
#> Predictions range: 0.0472 to 0.7025
#> Masked 116 of 225 grid points outside the data
land <- adaptive_landscape(prep, surf$model, c("total_length", "weight"),
                           simulation_n = 100, grid_n = 15)
#> Population mean grid ranges:
#>   total_length: -2.96 to 2.94
#>   weight: -3.12 to 4.85
#> Estimated within-population variance-covariance:
#>              total_length weight
#> total_length       1.0000 0.5839
#> weight             0.5839 1.0000
#> Calculating mean fitness for 225 grid points
#> Optimal population mean phenotype:
#>   total_length    weight
#> 7   -0.4336396 -3.121579
#> Mean fitness at optimum: 0.6725
#> The highest mean fitness is on the edge of the grid; the landscape may keep rising beyond it
#> 76% of simulated individuals fell outside the data (97% at the optimum)
#> The optimum rests on extrapolation: more than 25% of the population simulated there lies outside the data
comp <- compare_fitness_surfaces_data(surf, land, c("total_length", "weight"))
#> Warning: collinear points, trying to add some jitter to avoid colinearities!
#> Warning: success: collinearities reduced through jitter
#> Individual optimum:
#>    total_length    weight   fitness
#> 21   -0.6141753 -1.576041 0.7024572
#> Population optimum:
#>   total_length    weight   fitness
#> 7   -0.4336396 -3.121579 0.6724595
#> Distance between optima: 1.556
#> Summary statistics:
#>              Surface Fitness_Range_Min Fitness_Range_Max Fitness_Mean
#> 1 Correlated Fitness        0.08614618         0.7024572    0.4690910
#> 2 Adaptive Landscape        0.02567243         0.6724595    0.3835101
#>   Fitness_SD N_Points
#> 1  0.1845610      225
#> 2  0.1927433      225
#> Correlation between surfaces: 0.979
plot_fitness_surfaces_comparison(comp)
#> $side_by_side
#> Warning: Removed 116 rows containing non-finite outside the scale range
#> (`stat_contour_filled()`).

#> 
#> $overlay
#> Warning: Removed 116 rows containing non-finite outside the scale range
#> (`stat_contour()`).

#> 
```
