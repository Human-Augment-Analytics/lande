# Analyze linear selection gradients (beta)

Estimates linear selection gradients by OLS on relative fitness. Binary
fitness gets its p-values from a logistic GLM and count fitness from a
Poisson GLM (negative binomial if overdispersed), since the OLS tests
are not valid for either.

## Usage

``` r
analyze_linear_selection(
  data,
  fitness_col,
  trait_cols,
  fitness_type,
  binary_response_col = NULL,
  group = NULL,
  se_type = c("ols", "hc3")
)
```

## Arguments

- data:

  A data frame containing fitness and trait measurements.

- fitness_col:

  A string specifying the response column for the OLS gradient model
  (relative fitness).

- trait_cols:

  A character vector of trait column names.

- fitness_type:

  A string indicating the fitness type: `"binary"`, `"continuous"`,
  `"count"`, or `"proportion"`.

- binary_response_col:

  Optional string naming the raw fitness column (0/1 for binary, counts
  for count fitness) used for the GLM that supplies p-values. If `NULL`,
  `fitness_col` is treated as the raw outcome and relativised
  internally.

- group:

  Optional grouping column. With two or more groups the GLM that
  supplies the p-values gets an intercept for each, to match the
  standardising within groups; the least-squares fit has no group term.

- se_type:

  Standard errors of the least-squares gradients: `"ols"` (the default)
  or `"hc3"`, heteroscedasticity-consistent; see
  [`selection_coefficients()`](https://human-augment-analytics.github.io/lande/dev/reference/selection_coefficients.md).

## Value

A list containing the fitted models, summaries, ANOVA tables, and VIFs.

## Examples

``` r
prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
fit <- analyze_linear_selection(prep, "survival", c("total_length", "weight"), "binary")
extract_linear_coefficients(c("total_length", "weight"), fit)
#>           Term   Type Beta_Coefficient Standard_Error    P_Value
#> 1 total_length Linear      -0.17811220        0.09748 0.07276336
#> 2       weight Linear      -0.09992466        0.09748 0.30287778
```
