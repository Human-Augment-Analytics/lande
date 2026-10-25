# Calculate the Correlated Fitness Surface

Fits a model for individual fitness based on multiple individual
phenotypes (w ~ z1 + z2 + interactions). Traits MUST be standardized
BEFORE calling this function.

## Usage

``` r
correlated_fitness_surface(
  data,
  fitness_col,
  trait_cols,
  grid_n = 60,
  method = "auto",
  scale_traits = FALSE,
  group = NULL,
  group_effect = TRUE,
  k = NULL,
  mask = TRUE,
  too_far = NULL,
  bs = c("tp", "cr", "ps"),
  smoothing = c("REML", "GCV.Cp", "ML"),
  clamp = TRUE,
  level = 0.95,
  by_group = FALSE,
  count_family = c("poisson", "quasipoisson", "nb")
)
```

## Arguments

- data:

  A data frame containing fitness and trait measurements.

- fitness_col:

  A string specifying the name of the fitness column.

- trait_cols:

  A character vector of exactly length 2 specifying the trait column
  names.

- grid_n:

  Integer specifying the resolution of the prediction grid. Default is
  60.

- method:

  A string specifying the modeling method: `"auto"`, `"gam"`, or
  `"tps"`.

- scale_traits:

  Deprecated. Logical. Set to `FALSE` to avoid double standardization.

- group:

  Optional string naming a grouping column, such as species or year. The
  result's `groups` table then gives each group's mean trait values and
  the highest point of the surface within that group's own convex hull,
  and the plot functions draw both.

- group_effect:

  Logical; with `TRUE` (the default) the GAM includes `group` as a fixed
  effect and predicts the surface at the reference level, as the
  gradient models do for year or site; rows with no group label are one
  more level. With `FALSE` one surface is fitted to everyone and the
  group is only used for the `groups` table and the overlay, which is
  what a community surface of several species needs (Beausoleil et al.
  2023).

- k:

  Basis dimension for the GAM smooth. The default `NULL` sets it from
  the data as `min(30, max(10, floor(sqrt(n1 * n2))))`, where `n1` and
  `n2` are the numbers of distinct values of each trait; it is always
  capped at one less than the number of observations. Pass an integer to
  override.

- mask:

  Logical; if `TRUE` (the default) grid points outside the convex hull
  of the observed trait pairs get `NA` fitness, so the surface is only
  drawn where there are data. The grid column `.inside` records which
  points were kept, `.fit_all` holds the unmasked predictions, and the
  result's `hull` is the polygon the plot functions use to cover the
  outside.

- too_far:

  Optional distance rule for blanking the grid, on top of `mask`. Grid
  points farther than this from the nearest individual get `NA` fitness,
  with distances measured after scaling the grid to the unit square, as
  in [`mgcv::vis.gam`](https://rdrr.io/pkg/mgcv/man/vis.gam.html).
  Beausoleil et al. (2023) used 0.15. The default `NULL` applies no
  distance rule. The distance to the nearest individual is kept in the
  grid column `.dist`.

- bs:

  Basis for the GAM smooth: `"tp"` (thin plate, the default), `"cr"` or
  `"ps"`.

- smoothing:

  How the GAM's smoothing parameter is chosen: `"REML"` (the default),
  `"GCV.Cp"` or `"ML"`.

- clamp:

  Logical; with `method = "tps"` and `TRUE` (the default), predictions
  are held inside the range of the fitness type, 0 to 1 for survival and
  at least 0 for counts, as
  [`adaptive_landscape()`](https://human-augment-analytics.github.io/lande/dev/reference/adaptive_landscape.md)
  does. The GAM respects the range through its link and is unaffected.

- level:

  Confidence level of the band around a GAM surface. Default is 0.95.

- by_group:

  Logical; if `TRUE` fit a separate surface to each level of `group` and
  return a named list of surfaces, each with its own shape, grid and
  hull. The default `FALSE` fits one surface, with the group as a fixed
  effect or only marked on it (see `group_effect`).

- count_family:

  Family for count fitness in the GAM: `"poisson"` (the default);
  `"quasipoisson"`, which estimates the dispersion, so the standard
  errors and peak comparisons allow for overdispersed counts; or `"nb"`,
  a negative binomial. The smoothing parameter is chosen again with the
  dispersion estimated, so the fitted surface changes too, usually
  becoming smoother, and can change a lot where there are few
  individuals. The result's `dispersion` is the Pearson dispersion of
  the fit, and a Poisson fit warns when it is above 1.5 (a rule of
  thumb).

## Value

A list containing the fitted model, grid predictions, and metadata. With
`method = "gam"` the grid also carries `.se`, the standard error of the
fitted fitness, and `.fit_lo` and `.fit_hi`, the band at `level`, worked
out on the scale of the link and blanked where the fit is; the
thin-plate spline gives none. `peaks` lists the local maxima of the kept
surface, highest first, with `interior` `TRUE` when every neighbouring
cell is kept and lower, and `FALSE` when the maximum sits against the
edge of the data, which may only be where the data end. `original_data`
holds the rows that were analysed. With a `group` the `groups` data
frame has one row per group with its size, mean traits and the highest
point of the surface within its own range, flagged `peak_interior` when
that point is an interior maximum of the surface and `peak_edge` when it
lies at the edge of the data; with both `FALSE` the surface keeps rising
past the group's range. For count fitness `count_family` and
`dispersion` record the family used and the Pearson dispersion of the
fit.

## Details

The family follows the fitness column: 0/1 fitness gets a binomial
family, non-negative whole numbers with more than two values (lifespan
in years, offspring) a Poisson family with a log link or the one set by
`count_family`, and anything else a Gaussian family. With
`method = "auto"` binary and count fitness use the GAM and continuous
fitness the thin-plate spline. By default the GAM uses a thin-plate
smooth with the smoothing parameter chosen by REML; `bs` and `smoothing`
are there to match another study's smoother and do not apply to
`method = "tps"`. With `"cr"` or `"ps"`, which are one-dimensional
bases, the two traits enter as a tensor product smooth. Predicting a
fitted surface over the full rectangle of the grid extrapolates into
corners that no individual occupies; `mask = TRUE` leaves those blank
rather than showing a fitted value there. The hull still fills gaps
between separate clusters of individuals, such as several species on one
surface; `too_far` blanks those too.

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
surf$grid[which.max(surf$grid$.fit), ]
#>    total_length    weight      .fit       .se   .fit_lo   .fit_hi  .fit_all
#> 72   -0.5207948 -1.590066 0.7030432 0.1353279 0.3992416 0.8940022 0.7030432
#>      .se_all .fit_lo_all .fit_hi_all .inside
#> 72 0.1353279   0.3992416   0.8940022    TRUE

# two lakes of pupfish on one surface: cells far from any fish blank, each
# lake's mean and local peak marked, the lake kept out of the model
pup <- rbind(crescent_pond_pupfish, little_lake_pupfish)
pup <- pup[pup$density == "H", ]
prep2 <- prepare_selection_data(pup, "survival", c("nose", "noseangle"))
surf2 <- correlated_fitness_surface(prep2, "survival", c("nose", "noseangle"), grid_n = 30,
                                    too_far = 0.15, group = "lake", group_effect = FALSE)
#> Data type: binary; method: gam; n = 1671; k = 30
#> Grouping variable: lake (2 groups); one surface for all groups, the group only marks means and peaks
#> GAM fitting with 1671 observations
#>   Trying formula: main
#> Success with formula: main
#> Predictions range: 0.0843 to 0.2206
#> Masked 390 of 900 grid points outside the data
surf2$groups
#>   group   n    mean_nose mean_noseangle peak_nose peak_noseangle  peak_fit
#> 1    CP 796 -0.005094234    -0.01971555    4.6037       2.438193 0.2049404
#> 2    LL 875  0.004634298     0.01793552    4.6037       1.455267 0.2013855
#>   peak_interior peak_edge
#> 1         FALSE      TRUE
#> 2         FALSE      TRUE
```
