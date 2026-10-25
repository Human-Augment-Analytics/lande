# Prepare data for selection analysis

Cleans, standardizes, and calculates relative fitness for trait and
fitness data. Standardizations and relative fitness calculations can be
performed within groups.

## Usage

``` r
prepare_selection_data(
  data,
  fitness_col,
  trait_cols,
  standardize = TRUE,
  group = NULL,
  add_relative = TRUE,
  na_action = c("warn", "drop", "none"),
  name_relative = "relative_fitness"
)
```

## Arguments

- data:

  A data frame containing fitness and trait measurements.

- fitness_col:

  A string specifying the name of the fitness column.

- trait_cols:

  A character vector of trait column names.

- standardize:

  Logical indicating whether to standardize traits to mean 0 and SD 1.
  Default is `TRUE`.

- group:

  Optional string specifying a grouping variable.

- add_relative:

  Logical indicating whether to add a relative fitness column. Default
  is `TRUE`.

- na_action:

  A string specifying how to handle missing values: `"warn"`, `"drop"`,
  or `"none"`.

- name_relative:

  A string specifying the name for the relative fitness column. Default
  is `"relative_fitness"`.

## Value

A modified data frame ready for selection analysis.

## Examples

``` r
prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
head(prep[, c("survival", "relative_fitness", "total_length", "weight")])
#>   survival relative_fitness total_length     weight
#> 1        1         1.888889   -1.5569728 -0.6948140
#> 2        1         1.888889    0.1280269  0.9320675
#> 3        1         1.888889   -1.2761395  0.9320675
#> 4        1         1.888889   -1.5569728 -0.8303874
#> 5        1         1.888889   -0.9953062 -0.9659609
#> 6        1         1.888889    0.4088602  0.6609206

# standardised and relativised within each sex
by_sex <- prepare_selection_data(bumpus, "survival", "total_length", group = "sex")
#> Standardising and computing relative fitness within groups: 'sex'
tapply(by_sex$total_length, by_sex$sex, mean)
#>        female          male 
#> -1.249452e-15  2.333524e-15 
```
