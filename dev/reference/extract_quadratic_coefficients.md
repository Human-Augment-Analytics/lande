# Extract quadratic selection coefficients

Extract quadratic selection coefficients

## Usage

``` r
extract_quadratic_coefficients(trait_cols, results)
```

## Arguments

- trait_cols:

  A character vector of trait column names.

- results:

  A model results object returned by
  [`analyze_nonlinear_selection()`](https://human-augment-analytics.github.io/lande/dev/reference/analyze_nonlinear_selection.md).

## Value

A data frame with quadratic selection coefficients (doubled estimates)
and statistics.

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
```
