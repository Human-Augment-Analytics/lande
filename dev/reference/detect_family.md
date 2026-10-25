# Automatically detect fitness type and GLM family

Heuristically determines whether a fitness vector represents binary,
count, proportion, or continuous data.

## Usage

``` r
detect_family(y)
```

## Arguments

- y:

  A numeric vector representing fitness values.

## Value

A list containing `type` (string), `family` (GLM family object), and
`note` (string).

## Examples

``` r
detect_family(c(0, 1, 0, 0, 1, 1))
#> Warning: Binary fitness detected, but sample size < 10 (may be unstable)
#> $type
#> [1] "binary"
#> 
#> $family
#> 
#> Family: binomial 
#> Link function: logit 
#> 
#> 
#> $note
#> [1] "Binary fitness (0/1) detected. Use logistic GLM for p-values."
#> 
detect_family(bumpus$survival)$type
#> [1] "binary"
detect_family(rpois(50, 2))$type
#> [1] "count"
```
