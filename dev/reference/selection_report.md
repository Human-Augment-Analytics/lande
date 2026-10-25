# Standardised selection analysis report

Runs a Lande-Arnold analysis and returns one table of selection
differentials (S) and linear (beta), quadratic (gamma) and correlational
(gamma_ij) gradients. All estimates use traits standardised to mean 0,
SD 1 and relative fitness, so tables from different studies can be
compared.

## Usage

``` r
selection_report(
  data,
  fitness_col,
  trait_cols,
  fitness_type = c("auto", "binary", "count", "continuous"),
  standardize = TRUE,
  group = NULL,
  use_relative_for_fit = TRUE,
  include_differentials = TRUE,
  digits = 4,
  se_type = c("ols", "hc3")
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

- include_differentials:

  Logical; if `TRUE` (default) selection differentials S are included
  alongside the gradients.

- digits:

  Integer number of digits used when printing. Default is 4.

- se_type:

  Standard errors of the gradients: `"ols"`, the usual least-squares
  ones (the default), or `"hc3"`, leave-one-out standard errors that
  allow for residual spread changing with the traits. Mitchell-Olds and
  Shaw (1987) suggested the jackknife for selection gradients when the
  residuals are not normal; HC3 (MacKinnon and White 1985) gets much the
  same in closed form, from how far the estimates move when each
  individual is left out. In the package's simulation it recovered most
  of the coverage beta loses when the model leaves out curvature, from
  88 to 93% with normal traits, 67 to 86% with log-normal ones and 77 to
  91% with heavy-tailed symmetric ones. It did nothing for the coverage
  lost to estimating a skewed trait's SD, and with survival and counts
  its intervals covered a little less, down to 91% (see
  [`check_selection_assumptions()`](https://human-augment-analytics.github.io/lande/dev/reference/check_selection_assumptions.md)).
  For continuous fitness the p-values follow the chosen errors; for
  survival and counts they come from the GLM with either `se_type`.

## Value

A data frame of class `"selection_report"` with columns `Term`, `Type`,
`Estimate`, `Std_Error`, and `P_Value`.

## Details

Differentials and gradients are computed on the same individuals (those
with complete fitness and trait values) and on the same trait and
fitness scales. S is the population covariance (divides by n) while the
OLS gradient on sd-standardised traits corresponds to the sample
covariance (n - 1), so for a single trait the Differential and Linear
rows differ by the factor (n - 1) / n.

## Examples

``` r
selection_report(bumpus, "survival", c("total_length", "weight", "humerus"))
#> Warning: Collinear traits (VIF above 5) may inflate the standard errors
#> Selection analysis (standardised traits, relative fitness)
#> Fitness type: binary 
#> p-values are from a logistic model on the same terms
#> 
#>                    Term          Type Estimate Std_Error P_Value Sig
#>            total_length  Differential  -0.2347        NA      NA    
#>                  weight  Differential  -0.2024        NA      NA    
#>                 humerus  Differential   0.1629        NA      NA    
#>            total_length        Linear  -0.3153    0.0922  0.0016  **
#>                  weight        Linear  -0.2868    0.0952  0.0057  **
#>                 humerus        Linear   0.4833    0.0898  0.0000 ***
#>           total_length²     Quadratic  -0.1386    0.2088  0.3936    
#>                 weight²     Quadratic  -0.0098    0.1728  0.7222    
#>                humerus²     Quadratic   0.0056    0.1681  0.7165    
#>   total_length × weight Correlational  -0.0334    0.1609  0.5449    
#>  total_length × humerus Correlational  -0.0317    0.1191  0.7570    
#>        weight × humerus Correlational   0.0219    0.1560  0.5677    
#> 
#> Signif: *** 0.001  ** 0.01  * 0.05  . 0.1
```
