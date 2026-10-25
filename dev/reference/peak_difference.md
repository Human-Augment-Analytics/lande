# Compare the fitted fitness at two points of a surface

Takes two points on a GAM surface, given as group names or trait values,
and returns the difference in fitted fitness on the scale of the link,
with its standard error. Two peaks are separate only if each is higher
than the pass between them, the lowest point on the highest route from
one to the other, which `valley = TRUE` compares.

## Usage

``` r
peak_difference(
  surface,
  from,
  to,
  valley = FALSE,
  route = c("pass", "line"),
  n_path = 50
)
```

## Arguments

- surface:

  Output of
  [`correlated_fitness_surface()`](https://human-augment-analytics.github.io/lande/dev/reference/correlated_fitness_surface.md)
  fitted with `method = "gam"`.

- from, to:

  A group name, standing for that group's highest point in
  `surface$groups`, or a numeric vector of the two trait values.

- valley:

  Logical; also find the valley between the two points and compare each
  end with it. Default is `FALSE`.

- route:

  How the valley is found: `"pass"` (the default), the pass over the
  kept cells of the surface's grid, or `"line"`, the lowest point on the
  straight line between the two points, which is how Beausoleil et
  al. (2023) measured valley depths. A straight line can cross a trough
  that a curving ridge goes round, so it can make two peaks look more
  separate than they are.

- n_path:

  Number of points along the straight line with `route = "line"`.
  Default is 50.

## Value

A data frame with one row per comparison: the fitted fitness at the two
points (`fit_a`, `fit_b`), their difference on the link scale, its
standard error, `z` and the two-sided `p_value` (`NA` for the valley
rows). With `valley = TRUE` the attribute `"valley"` holds the trait
values of the pass or of the lowest point on the line.

## Details

Differences are on the scale of the link, log fitness for counts and log
odds for survival. The standard error comes from the covariance of the
model's coefficients, corrected for smoothing parameter uncertainty when
the fit was by REML or ML. The points are taken as fixed, although the
peaks and the valley were found on the same fitted surface, so the
comparisons describe the fitted surface and are not planned tests. The
valley is the lowest point on the route, so its differences are biased
upwards and their z values too large, which is why the valley rows have
no p-value. The pass is found cell by cell over the grid, joining each
kept cell to its eight neighbours, from the two kept cells nearest the
points in grid steps; its precision follows `grid_n`, and a point far
from any kept cell gets a message. If the route never drops below the
cell a point starts from, or the lower point is no higher than the pass,
there is no dip to compare. The straight line ignores the mask. For
overdispersed counts fit the surface with
`count_family = "quasipoisson"` first, since Poisson standard errors are
then too small.

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
peak_difference(surf, c(-1, -1), c(1, 1))
#>          comparison     fit_a     fit_b difference        se        z
#> 1 (-1, -1) - (1, 1) 0.6883016 0.3811798   1.276735 0.6329478 2.017125
#>      p_value
#> 1 0.04368246
```
