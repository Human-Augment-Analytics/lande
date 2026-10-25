# Extract linear selection coefficients

Extract linear selection coefficients

## Usage

``` r
extract_linear_coefficients(trait_cols, results)
```

## Arguments

- trait_cols:

  A character vector of trait column names.

- results:

  A model results object returned by
  [`analyze_linear_selection()`](https://human-augment-analytics.github.io/lande/dev/reference/analyze_linear_selection.md).

## Value

A data frame with linear selection coefficients and statistics.

## Examples

``` r
prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
fit <- analyze_linear_selection(prep, "survival", c("total_length", "weight"), "binary")
extract_linear_coefficients(c("total_length", "weight"), fit)
#>           Term   Type Beta_Coefficient Standard_Error    P_Value
#> 1 total_length Linear      -0.17811220        0.09748 0.07276336
#> 2       weight Linear      -0.09992466        0.09748 0.30287778
```
