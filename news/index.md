# Changelog

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
  [`univariate_spline()`](https://human-augment-analytics.github.io/lande/reference/univariate_spline.md),
  [`correlated_fitness_surface()`](https://human-augment-analytics.github.io/lande/reference/correlated_fitness_surface.md)
  and
  [`temporal_landscape()`](https://human-augment-analytics.github.io/lande/reference/temporal_landscape.md)
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
  [`run_app()`](https://human-augment-analytics.github.io/lande/reference/run_app.md),
  which writes out the R calls that repeat an analysis with the current
  settings.
- [`canonical_analysis()`](https://human-augment-analytics.github.io/lande/reference/canonical_analysis.md)
  rotates the gamma matrix to its canonical axes, with double-regression
  standard errors, permutation p-values (Reynolds et al.
  2010. and an optional bootstrap, and
        [`plot_canonical_axes()`](https://human-augment-analytics.github.io/lande/reference/plot_canonical_axes.md)
        draws fitness along them.
- Standard errors and a band on the GAM fitness surface, drawn with the
  plots’ `uncertainty` option, and
  [`peak_difference()`](https://human-augment-analytics.github.io/lande/reference/peak_difference.md)
  to compare two points of a surface, or each with the pass between them
  (`route = "line"` for the straight line).
- The adaptive landscape reports how much of each simulated population
  falls outside the data, and flags an optimum on the edge of the grid;
  the plot can mark where the population leaves the data.
- A `clamp` option on the landscape and the thin-plate surface holds
  predicted fitness inside the range of the fitness type; on by default.
