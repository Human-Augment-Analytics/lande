# Extract correlational selection coefficients

Extract correlational selection coefficients

## Usage

``` r
extract_interaction_coefficients(trait_cols, results)
```

## Arguments

- trait_cols:

  A character vector of trait column names.

- results:

  A model results object returned by
  [`analyze_nonlinear_selection()`](https://human-augment-analytics.github.io/lande/dev/reference/analyze_nonlinear_selection.md).

## Value

A data frame with correlational (interaction) selection coefficients and
statistics.

## Examples

``` r
prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
fit <- analyze_nonlinear_selection(prep, "survival", c("total_length", "weight"), "binary")
#> there are higher-order terms (interactions) in this model
#> consider setting type = 'predictor'; see ?vif
#> Warning: Collinear traits (VIF above 5) may inflate the standard errors
extract_interaction_coefficients(c("total_length", "weight"), fit)
#>                    Term          Type Beta_Coefficient Standard_Error   P_Value
#> 1 total_length × weight Correlational      -0.06949785      0.1585041 0.5723638
```
