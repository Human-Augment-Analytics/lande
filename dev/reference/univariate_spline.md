# Estimate univariate correlated fitness function

This function calculates a univariate correlated fitness function using
a GAM smooth term. The formula is `w ~ f(z)`. Traits should be
standardized BEFORE calling this function.

## Usage

``` r
univariate_spline(
  data,
  fitness_col,
  trait_col,
  fitness_type = c("auto", "binary", "count", "continuous"),
  group = NULL,
  relative_col = NULL,
  k = 10,
  bs = c("cr", "tp", "ps"),
  smoothing = c("GCV.Cp", "REML", "ML"),
  bootstrap = FALSE,
  n_boot = 1000,
  by_group = FALSE,
  count_family = c("poisson", "quasipoisson", "nb"),
  cl = NULL
)
```

## Arguments

- data:

  A data frame containing fitness and trait measurements.

- fitness_col:

  A string specifying the name of the fitness column.

- trait_col:

  A string specifying the name of the trait column (must be numeric and
  standardized).

- fitness_type:

  A string indicating the fitness type: `"auto"` (detect from the data,
  the default), `"binary"`, `"count"`, or `"continuous"`. Binary and
  count fitness are fitted on the raw values with a binomial or Poisson
  family; continuous fitness on relative fitness with a Gaussian one.

- group:

  Optional string specifying a grouping variable. If provided, group
  fixed effects are included (none when it has only one level); rows
  with no group label are one more level.

- relative_col:

  Optional string naming a pre-computed relative fitness column to use
  for continuous fitness (e.g. one produced within groups by
  `prepare_selection_data`).

- k:

  Integer specifying the basis dimension for the smooth term. Default is
  10.

- bs:

  Spline basis: `"cr"` (cubic regression spline, the default), `"tp"`
  (thin plate) or `"ps"` (P-spline).

- smoothing:

  How the smoothing parameter is chosen: `"GCV.Cp"` (generalised
  cross-validation, the default, which mgcv applies in its UBRE form to
  binary and Poisson fitness and as plain GCV to quasi-Poisson; a
  negative binomial fit uses REML in its place), `"REML"` or `"ML"`. The
  result's `spline_type` names the criterion mgcv used.

- bootstrap:

  Logical; if `TRUE` the 95% ribbon is obtained by resampling
  individuals and refitting (Schluter 1988). The default `FALSE` uses
  the parametric Wald interval, which is fast; the bootstrap refits the
  spline `n_boot` times.

- n_boot:

  Integer number of bootstrap resamples used when `bootstrap = TRUE`.
  Default is 1000.

- by_group:

  Logical; if `TRUE` fit each level of `group` on its own rows and
  return a named list of fits, one per group. The default `FALSE` fits
  one curve with the group as a fixed effect. Prepare the data with the
  same `group` first, so that each group is standardised on its own.

- count_family:

  Family for count fitness: `"poisson"` (the default), `"quasipoisson"`
  or `"nb"`, as in
  [`correlated_fitness_surface()`](https://human-augment-analytics.github.io/lande/dev/reference/correlated_fitness_surface.md);
  the smoothing parameter is chosen under whichever is used. The
  result's `dispersion` is the Pearson dispersion of the fit, and a
  Poisson fit warns when it is above 1.5.

- cl:

  Optional cluster from
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html)
  to run the bootstrap refits on. The default `NULL` runs them one after
  another in this session. With a cluster the resamples are still drawn
  here, in the same order, and only the refits go to the workers, which
  fit under this session's random number generator kinds and its options
  for contrasts, matrix products and warnings. Their warnings and
  messages are raised here in the order the serial loop raises them, so
  the ribbon, the warnings, `options(warn = 2)` and the state of the
  random number generator afterwards (or when a handler stops the call
  at a warning) are the same as without a cluster. A negative binomial
  fit always runs here, as each refit starts from the theta of the one
  before. The caller starts and stops the cluster. The workers fit with
  mgcv, so they need the same version of it as this session on their
  library paths (for example
  `parallel::clusterCall(cl, function(p) .libPaths(p), .libPaths())`);
  lande itself need not be loaded on them.

## Value

A list of class `"univariate_fitness"` containing the fitted GAM model,
a prediction grid, and metadata.

## Details

By default the fitness function is a penalised cubic regression spline
with the smoothing parameter chosen by generalised cross-validation
(UBRE for binary and Poisson fitness), following Schluter (1988). `bs`
and `smoothing` set the family of curves and the criterion, to match the
smoother of another study; the amount of smoothing is still estimated
from the data, though the criterion can change it a good deal. If mgcv's
check suggests the basis dimension was too small a warning says so; the
check works on the residuals, which for survival carry too little to
test, so it is skipped there. With `bootstrap = TRUE` the result depends
on the random seed; call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first for a
reproducible ribbon.

## Examples

``` r
prep <- prepare_selection_data(bumpus, "survival", "total_length")
uni <- univariate_spline(prep, "survival", "total_length")
#> Fitness type detected: binary
plot_univariate_fitness(uni, "total_length")


# a thin-plate basis with REML, to match another study's smoother
univariate_spline(prep, "survival", "total_length", bs = "tp", smoothing = "REML")$spline_type
#> Fitness type detected: binary
#> [1] "thin-plate spline (REML)"
```
