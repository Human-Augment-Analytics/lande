# Bootstrap selection gradients

Resamples individuals with replacement and re-estimates the selection
gradients on each resample to obtain bootstrap standard errors and
percentile confidence intervals for every coefficient (beta, gamma,
gamma_ij).

## Usage

``` r
bootstrap_selection(
  data,
  fitness_col,
  trait_cols,
  fitness_type = c("auto", "binary", "count", "continuous"),
  standardize = TRUE,
  group = NULL,
  use_relative_for_fit = TRUE,
  n_boot = 1000,
  conf = 0.95
)
```

## Arguments

- data:

  A data frame containing fitness and trait measurements.

- fitness_col:

  A string specifying the name of the fitness column.

- trait_cols:

  A character vector of trait column names.

- fitness_type:

  A string indicating the fitness type: `"auto"`, `"binary"`, `"count"`,
  or `"continuous"`. Binary and count fitness take their p-values from a
  GLM (logistic, or Poisson and negative binomial) on the raw values;
  the gradients always come from OLS on relative fitness.

- standardize:

  Logical indicating whether to standardize traits to mean 0 and SD 1.
  Default is `TRUE`.

- group:

  Optional string specifying a grouping variable (e.g., "year", "site").
  Traits and relative fitness are standardised within each group, and
  for binary and count fitness the GLM that supplies the p-values gets a
  separate intercept for each group.

- use_relative_for_fit:

  Logical; if `TRUE` (default) the gradients are estimated on relative
  fitness \\W / \bar{W}\\ for every fitness type, as in Lande & Arnold
  (1983). Set `FALSE` only to reproduce coefficients on the absolute
  fitness scale.

- n_boot:

  Integer number of bootstrap resamples. Default is 1000.

- conf:

  Confidence level for the percentile interval. Default is 0.95.

## Value

A data frame with one row per coefficient and columns `Term`, `Type`,
`Estimate` (point estimate on the full data), `Boot_SE`, `CI_lower`,
`CI_upper`, `P_Value`, and `N_Boot` (usable resamples for that
coefficient). `P_Value` is the parametric p-value from the
point-estimate fit (the OLS t-test for continuous fitness, the Wald test
of the logistic GLM for binary fitness or of the Poisson or negative
binomial GLM for counts), not a bootstrap p-value.

## Details

Each resample goes through the whole procedure again: with
`standardize = TRUE` the traits are restandardised, and fitness is made
relative to the resample's own mean, before both least-squares models
are refitted, so the intervals carry the uncertainty in the
standardisation as well as in the fit. Only those two models are
refitted on each resample; the GLMs that supply the p-values (logistic
for binary fitness, Poisson or negative binomial for counts) are fitted
once, on the full data, so a GLM that would fail on some resample cannot
remove it. Resamples whose least-squares fit fails are dropped, and
`N_Boot` counts the ones that were used. If the column names keep this
shortcut from reproducing the point estimate exactly (a trait named like
an R expression, say), each resample is refitted with all of
[`selection_coefficients()`](https://human-augment-analytics.github.io/lande/dev/reference/selection_coefficients.md)
instead, which is slower and drops a resample whose GLM fails.

When `group` is given, individuals are resampled within each group so
every resample keeps the original group sizes. A resample in which a
trait has no variance within a group (so it could only be centred) is
discarded rather than fitted on a degenerate value. The point estimate
is the ordinary
[`selection_coefficients()`](https://human-augment-analytics.github.io/lande/dev/reference/selection_coefficients.md)
fit on the full data and its warnings are reported as usual.

With continuous fitness the intervals for beta can be too narrow when
the traits have heavy tails (skewed or not) and the linear model leaves
out curvature, though less so than the parametric ones (see
[`check_selection_assumptions()`](https://human-augment-analytics.github.io/lande/dev/reference/check_selection_assumptions.md)).
For survival and counts in the package's simulation they covered a
little less than the parametric ones, 92 to 93% for beta. Refitting on a
transformed trait changes the scale of the gradients and is only a
check.

## Examples

``` r
set.seed(1)
bootstrap_selection(bumpus, "survival", c("total_length", "weight"), n_boot = 50)
#> Warning: Collinear traits (VIF above 5) may inflate the standard errors
#>                    Term          Type    Estimate    Boot_SE   CI_lower
#> 1          total_length        Linear -0.17811220 0.10303749 -0.4021452
#> 2                weight        Linear -0.09992466 0.08007062 -0.2573204
#> 3         total_length²     Quadratic -0.23342057 0.26857154 -0.7271855
#> 4               weight²     Quadratic  0.02096490 0.17479658 -0.2530171
#> 5 total_length × weight Correlational -0.06949785 0.16260489 -0.3050220
#>      CI_upper    P_Value N_Boot
#> 1 0.005323691 0.07276336     50
#> 2 0.032591335 0.30287778     50
#> 3 0.174355313 0.24402273     50
#> 4 0.329534214 0.93873420     50
#> 5 0.271905000 0.57236384     50
```
