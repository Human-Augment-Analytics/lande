# Changelog

## lande (development version)

- On a server, the app limits a run to 10 traits, 5000 individuals and
  20 groups, and draws fewer bootstrap resamples for a large fit, so one
  Run can’t hold up everyone else.
  [`run_app()`](https://human-augment-analytics.github.io/lande/dev/reference/run_app.md)
  has no limits.
- [`plot_correlated_fitness_3d()`](https://human-augment-analytics.github.io/lande/dev/reference/plot_correlated_fitness_3d.md)
  draws the fitted fitness surface in 3D as a mesh of grid squares
  outlined in dark grey, only where there are data, with its interior
  peaks marked and light shading that `shade` sets.
  `color_palette = "jet"` gives the rainbow ramp of many published
  surfaces; viridis stays the default as it reads better for
  colour-blind viewers. The app’s fitness surface tab can show it too,
  with the group means and high points, and a jet colour theme for the
  2D and 3D views.
- The app is much faster. Run shows the gradients at once, each fit is
  kept so a change refits only what it affects, and the spline band runs
  on worker processes (`options(lande.workers = 0)` turns them off). The
  landscape now starts on the lattice method.
- [`univariate_spline()`](https://human-augment-analytics.github.io/lande/dev/reference/univariate_spline.md)
  gains `cl` to run the bootstrap refits on a cluster, with the same
  results.
- [`bootstrap_selection()`](https://human-augment-analytics.github.io/lande/dev/reference/bootstrap_selection.md)
  refits only the least-squares gradients on each resample and runs 40
  to 55 times faster. A resample is no longer dropped when its p-value
  model would fail to fit.
- [`selection_coefficients()`](https://human-augment-analytics.github.io/lande/dev/reference/selection_coefficients.md)
  no longer fails when a least-squares fit is exact.
- [`adaptive_landscape()`](https://human-augment-analytics.github.io/lande/dev/reference/adaptive_landscape.md)
  gains `method = "lattice"`, the same draws as before but 4 to 10 times
  faster, and `method = "convolve"`, a Gaussian convolution with no
  simulation noise. The default is unchanged.
- The app reads semicolon, tab and UTF-16 CSV files and Excel’s Windows
  encoding, and says why a file can’t be read.
- The app can download the whole analysis as one zip, with the data and
  an R script that reruns it.
- The app has a new layout built around its results, each with its own
  downloads.

## lande 0.1.1

- [`adaptive_landscape()`](https://human-augment-analytics.github.io/lande/dev/reference/adaptive_landscape.md)
  holds a numeric group, such as year, at the same reference level as
  [`correlated_fitness_surface()`](https://human-augment-analytics.github.io/lande/dev/reference/correlated_fitness_surface.md).
  It took the median level of the raw column, while the surface took the
  most common level, so a landscape from a surface with `group = "year"`
  could be drawn at a different year from the surface.
- The app no longer needs plotly to start: without it, the 3D landscape
  falls back to the fields plot as intended.
- The finch_community help page cites the authors’ CC BY 4.0 data
  deposit on Borealis as its source.
- The title is now “Standardised Analyses of Natural Selection Data”.
- The minimum versions of dplyr (1.1.0) and testthat (3.1.7) are set.
- Smaller fixes to the help pages, the examples and the vignette.

## lande 0.1.0

CRAN release: 2026-10-16

First release.

- Selection differentials and linear, quadratic and correlational
  gradients on relative fitness, with p-values from the logistic model
  for survival and from a Poisson or negative binomial model for counts.
- Cubic spline fitness functions with a bootstrap band, two-trait
  fitness surfaces drawn only where there are data, and adaptive
  landscapes of mean fitness for one or two traits.
- `count_family = "quasipoisson"` or `"nb"` in
  [`univariate_spline()`](https://human-augment-analytics.github.io/lande/dev/reference/univariate_spline.md),
  [`correlated_fitness_surface()`](https://human-augment-analytics.github.io/lande/dev/reference/correlated_fitness_surface.md)
  and
  [`temporal_landscape()`](https://human-augment-analytics.github.io/lande/dev/reference/temporal_landscape.md)
  for overdispersed counts; the spline and the surface report the
  dispersion.
- `by_group = TRUE` in the spline and the surface fits each group on its
  own and returns a list of fits.
- Several groups on one surface, with each group’s mean and highest
  point marked, and interior peaks told from maxima at the edge of the
  data.
- Fitness functions and landscapes fitted period by period.
- Bootstrap standard errors and percentile intervals for the gradients,
  and leave-one-out (HC3) standard errors with `se_type = "hc3"`.
- Assumption checks in one table: multivariate normality of the traits,
  VIF, rows per term, residuals or dispersion of the gradient model.
- Bumpus’s sparrows, Martin’s pupfish and two sets of Galapagos finches
  bundled as datasets, with CSV copies of the pupfish and finch data in
  `extdata`.
- A Shiny app that runs the analyses in a browser, opened with
  [`run_app()`](https://human-augment-analytics.github.io/lande/dev/reference/run_app.md),
  which writes out the R calls that repeat an analysis with the current
  settings.
- [`canonical_analysis()`](https://human-augment-analytics.github.io/lande/dev/reference/canonical_analysis.md)
  rotates the gamma matrix to its canonical axes, with double-regression
  standard errors, permutation p-values (Reynolds et al.
  2010. and an optional bootstrap, and
        [`plot_canonical_axes()`](https://human-augment-analytics.github.io/lande/dev/reference/plot_canonical_axes.md)
        draws fitness along them.
- Standard errors and a band on the GAM fitness surface, drawn with the
  plots’ `uncertainty` option, and
  [`peak_difference()`](https://human-augment-analytics.github.io/lande/dev/reference/peak_difference.md)
  to compare two points of a surface, or each with the pass between them
  (`route = "line"` for the straight line).
- The adaptive landscape reports how much of each simulated population
  falls outside the data, and flags an optimum on the edge of the grid;
  the plot can mark where the population leaves the data.
- A `clamp` option on the landscape and the thin-plate surface holds
  predicted fitness inside the range of the fitness type; on by default.
