# Calculate selection differential (S)

Computes the selection differential S = Cov(z, w) for a single trait.

## Usage

``` r
selection_differential(
  data,
  fitness_col,
  trait_col,
  standardized = TRUE,
  use_relative = TRUE,
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

  A string specifying the name of the trait column.

- standardized:

  Logical indicating whether the trait is already standardized. If
  `FALSE`, it will be standardized.

- use_relative:

  Logical indicating whether to use relative fitness.

- group:

  Optional grouping variable string. Calculates S within each group if
  provided.

- return_grouped:

  Logical indicating whether to return a data frame with S per group
  (`TRUE`) or overall mean S (`FALSE`).

## Value

A single numeric value representing overall S, or a data frame with S
per group if `return_grouped = TRUE`.

## Examples

``` r
prep <- prepare_selection_data(bumpus, "survival", "total_length")
selection_differential(prep, "survival", "total_length")
#> [1] -0.2347161
selection_differential(prep, "survival", "total_length", group = "sex", return_grouped = TRUE)
#> # A tibble: 2 × 3
#>   sex         S     n
#>   <fct>   <dbl> <int>
#> 1 female -0.168    49
#> 2 male   -0.329    87
```
