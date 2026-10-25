# Disruptive or stabilising selection on one trait

Estimates the linear (beta) and quadratic (gamma) selection gradients of
a single trait. Gamma is positive under disruptive selection and
negative under stabilising selection.

## Usage

``` r
analyze_disruptive_selection(
  data,
  fitness_col,
  trait_col,
  fitness_type = c("auto", "binary", "count", "continuous"),
  standardize = TRUE,
  group = NULL,
  return_grouped = FALSE
)
```

## Arguments

- data:

  A data frame containing fitness and trait measurements.

- fitness_col:

  A string specifying the name of the fitness column.

- trait_col:

  A string specifying the name of the single trait column.

- fitness_type:

  A string indicating the fitness type: `"auto"` (detect from the data),
  `"binary"`, `"count"` or `"continuous"`. Default is `"auto"`.

- standardize:

  Logical indicating whether to standardize the trait to mean 0 and
  SD 1. Default is `TRUE`.

- group:

  Optional string specifying a grouping variable; standardisation and
  relative fitness are then computed within each group.

- return_grouped:

  Logical; if `TRUE` and `group` is given, the gradients are estimated
  separately for each group and returned with a `Group` column. Default
  is `FALSE`.

## Value

A data frame with one row per gradient (`Term`, `Type`,
`Beta_Coefficient`, `Standard_Error`, `P_Value`, `Variance`).

## Examples

``` r
analyze_disruptive_selection(bumpus, "survival", "total_length")
#>            Term      Type Beta_Coefficient Standard_Error     P_Value
#> 1  total_length    Linear       -0.2364547     0.07915422 0.004487591
#> 2 total_length² Quadratic       -0.3245696     0.13732470 0.019694542
#>      Variance
#> 1 0.006265391
#> 2 0.018858072
analyze_disruptive_selection(bumpus, "survival", "total_length",
                             group = "sex", return_grouped = TRUE)
#>                   Term      Type Beta_Coefficient Standard_Error      P_Value
#> male.1    total_length    Linear       -0.3694482     0.08243345 0.0002002952
#> male.2   total_length² Quadratic       -0.1042183     0.13610059 0.7263041851
#> female.1  total_length    Linear       -0.1672317     0.16841847 0.3189315357
#> female.2 total_length² Quadratic       -0.3448404     0.36946919 0.3292407999
#>             Variance  Group
#> male.1   0.006795274   male
#> male.2   0.018523371   male
#> female.1 0.028364780 female
#> female.2 0.136507482 female
```
