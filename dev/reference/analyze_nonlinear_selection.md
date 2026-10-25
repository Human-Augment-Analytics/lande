# Analyze nonlinear selection gradients (gamma)

Estimates quadratic and correlational selection gradients by OLS on
relative fitness.

## Usage

``` r
analyze_nonlinear_selection(
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

  Optional grouping column, used as in
  [`analyze_linear_selection()`](https://human-augment-analytics.github.io/lande/dev/reference/analyze_linear_selection.md):
  a separate intercept for each group in the GLM that supplies the
  p-values.

- se_type:

  Standard errors of the least-squares gradients: `"ols"` (the default)
  or `"hc3"`, heteroscedasticity-consistent; see
  [`selection_coefficients()`](https://human-augment-analytics.github.io/lande/dev/reference/selection_coefficients.md).

## Value

A list containing the fitted nonlinear models, summaries, ANOVA tables,
and VIFs.

## Examples

``` r
prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
fit <- analyze_nonlinear_selection(prep, "survival", c("total_length", "weight"), "binary")
#> there are higher-order terms (interactions) in this model
#> consider setting type = 'predictor'; see ?vif
#> Warning: Collinear traits (VIF above 5) may inflate the standard errors
extract_quadratic_coefficients(c("total_length", "weight"), fit)
#>            Term      Type Beta_Coefficient Standard_Error   P_Value
#> 1 total_length² Quadratic       -0.2334206      0.2111590 0.2440227
#> 2       weight² Quadratic        0.0209649      0.1628337 0.9387342
extract_interaction_coefficients(c("total_length", "weight"), fit)
#>                    Term          Type Beta_Coefficient Standard_Error   P_Value
#> 1 total_length × weight Correlational      -0.06949785      0.1585041 0.5723638
```
