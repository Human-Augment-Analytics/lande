# Calculate selection coefficients

Main wrapper function for selection analysis following Lande & Arnold
(1983). Extracts linear (beta), quadratic (gamma), and correlational
selection gradients.

## Usage

``` r
selection_coefficients(
  data,
  fitness_col,
  trait_cols,
  fitness_type = c("auto", "binary", "count", "continuous"),
  standardize = TRUE,
  group = NULL,
  use_relative_for_fit = TRUE,
  return_grouped = FALSE,
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

- return_grouped:

  Logical indicating whether to return results grouped if a `group` is
  specified.

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

A data frame containing selection coefficients (Term, Type,
Beta_Coefficient, Standard_Error, P_Value, Variance), with the standard
errors used in the attribute `"se_type"`.

## References

MacKinnon, J. G. and White, H. (1985) Some heteroskedasticity-consistent
covariance matrix estimators with improved finite sample properties.
Journal of Econometrics 29, 305-325. Mitchell-Olds, T. and Shaw, R. G.
(1987) Regression analysis of natural selection: statistical inference
and biological interpretation. Evolution 41, 1149-1161.

## Examples

``` r
selection_coefficients(bumpus, "survival", c("total_length", "weight"), fitness_type = "binary")
#> there are higher-order terms (interactions) in this model
#> consider setting type = 'predictor'; see ?vif
#> Warning: Collinear traits (VIF above 5) may inflate the standard errors
#>                    Term          Type Beta_Coefficient Standard_Error
#> 1          total_length        Linear      -0.17811220      0.0974800
#> 2                weight        Linear      -0.09992466      0.0974800
#> 3         total_length²     Quadratic      -0.23342057      0.2111590
#> 4               weight²     Quadratic       0.02096490      0.1628337
#> 5 total_length × weight Correlational      -0.06949785      0.1585041
#>      P_Value    Variance
#> 1 0.07276336 0.009502351
#> 2 0.30287778 0.009502351
#> 3 0.24402273 0.044588139
#> 4 0.93873420 0.026514809
#> 5 0.57236384 0.025123549

# one set of gradients per sex, each sex standardised on its own
selection_coefficients(bumpus, "survival", c("total_length", "weight"),
                       group = "sex", return_grouped = TRUE)
#> there are higher-order terms (interactions) in this model
#> consider setting type = 'predictor'; see ?vif
#> there are higher-order terms (interactions) in this model
#> consider setting type = 'predictor'; see ?vif
#>                           Term          Type Beta_Coefficient Standard_Error
#> male.1            total_length        Linear    -0.3407483627     0.09731904
#> male.2                  weight        Linear    -0.0545617939     0.09731904
#> male.3           total_length²     Quadratic    -0.0273816426     0.19241802
#> male.4                 weight²     Quadratic     0.0462210677     0.17271044
#> male.5   total_length × weight Correlational    -0.0762079316     0.15832913
#> female.1          total_length        Linear     0.0004585805     0.20517997
#> female.2                weight        Linear    -0.2876716237     0.20517997
#> female.3         total_length²     Quadratic    -0.2985821975     0.48947090
#> female.4               weight²     Quadratic     0.0333986785     0.32247945
#> female.5 total_length × weight Correlational    -0.0520198090     0.33611268
#>              P_Value    Variance  Group
#> male.1   0.001978417 0.009470995   male
#> male.2   0.606430455 0.009470995   male
#> male.3   0.832288575 0.037024693   male
#> male.4   0.751956548 0.029828895   male
#> male.5   0.544449522 0.025068112   male
#> female.1 0.962802935 0.042098822 female
#> female.2 0.159163810 0.042098822 female
#> female.3 0.478153622 0.239581761 female
#> female.4 0.599203577 0.103992995 female
#> female.5 0.968497077 0.112971735 female
```
