# lande 0.1.1

* `adaptive_landscape()` holds a numeric group, such as year, at the same
  reference level as `correlated_fitness_surface()`. It took the median level
  of the raw column, while the surface took the most common level, so a
  landscape from a surface with `group = "year"` could be drawn at a different
  year from the surface.
* The app no longer needs plotly to start: without it, the 3D landscape falls
  back to the fields plot as intended.
* The finch_community help page cites the authors' CC BY 4.0 data deposit on
  Borealis as its source.
* The title is now "Standardised Analyses of Natural Selection Data".
* The minimum versions of dplyr (1.1.0) and testthat (3.1.7) are set.
* Smaller fixes to the help pages, the examples and the vignette.

# lande 0.1.0

First release.

* Selection differentials and linear, quadratic and correlational gradients on
  relative fitness, with p-values from the logistic model for survival and from
  a Poisson or negative binomial model for counts.
* Cubic spline fitness functions with a bootstrap band, two-trait fitness
  surfaces drawn only where there are data, and adaptive landscapes of mean
  fitness for one or two traits.
* `count_family = "quasipoisson"` or `"nb"` in `univariate_spline()`,
  `correlated_fitness_surface()` and `temporal_landscape()` for overdispersed
  counts; the spline and the surface report the dispersion.
* `by_group = TRUE` in the spline and the surface fits each group on its own
  and returns a list of fits.
* Several groups on one surface, with each group's mean and highest point
  marked, and interior peaks told from maxima at the edge of the data.
* Fitness functions and landscapes fitted period by period.
* Bootstrap standard errors and percentile intervals for the gradients, and
  leave-one-out (HC3) standard errors with `se_type = "hc3"`.
* Assumption checks in one table: multivariate normality of the traits, VIF,
  rows per term, residuals or dispersion of the gradient model.
* Bumpus's sparrows, Martin's pupfish and two sets of Galapagos finches
  bundled as datasets, with CSV copies of the pupfish and finch data in
  `extdata`.
* A Shiny app that runs the analyses in a browser, opened with `run_app()`,
  which writes out the R calls that repeat an analysis with the current
  settings.
* `canonical_analysis()` rotates the gamma matrix to its canonical axes, with
  double-regression standard errors, permutation p-values (Reynolds et al.
  2010) and an optional bootstrap, and `plot_canonical_axes()` draws fitness
  along them.
* Standard errors and a band on the GAM fitness surface, drawn with the plots'
  `uncertainty` option, and `peak_difference()` to compare two points of a
  surface, or each with the pass between them (`route = "line"` for the
  straight line).
* The adaptive landscape reports how much of each simulated population falls
  outside the data, and flags an optimum on the edge of the grid; the plot can
  mark where the population leaves the data.
* A `clamp` option on the landscape and the thin-plate surface holds predicted
  fitness inside the range of the fitness type; on by default.
