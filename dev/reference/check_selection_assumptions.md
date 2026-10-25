# Check the assumptions behind a selection analysis

Runs the checks a Lande and Arnold analysis rests on and returns them in
one table: multivariate normality of the traits, normality of each
trait, collinearity, rows per model term, and for the gradient models a
residual normality and heteroscedasticity test (continuous fitness), a
separation check (binary fitness) or the dispersion ratio (count
fitness). Nothing is enforced; the table is there to report alongside
the gradients, as the protocol of Palacio et al. (2019) asks.

## Usage

``` r
check_selection_assumptions(
  data,
  fitness_col,
  trait_cols,
  fitness_type = c("auto", "binary", "count", "continuous"),
  standardize = TRUE,
  group = NULL
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

## Value

A data frame of class `"selection_assumptions"` with one row per check:
`check`, `statistic`, `p_value` and `note`.

## Details

The regression gradients equal the covariance-adjusted differentials,
\\P^{-1} S\\, by least squares whatever the trait distribution. Normal
traits let them be read as the average slope and curvature of the
fitness surface (Lande and Arnold 1983; Morrissey and Sakrejda 2013) and
tie gamma to the change in the phenotypic covariance. Fitness doesn't
need to be normal, and normality doesn't change how the gradients are
computed. Mardia's (1970) skewness and kurtosis test multivariate
normality, and each trait is also tested on its own with Shapiro and
Wilk's test when there are 5000 rows or fewer.

In the package's simulation (validation/simulation.R in the GitHub
repository), on a curved fitness surface, nominal 95% intervals for beta
covered 88% of the time with normal traits, 67% with log-normal ones and
77% with symmetric heavy-tailed ones (92, 80 and 91% by bootstrap). Most
of that loss comes from the curvature, and HC3 standard errors
(`selection_coefficients(se_type = "hc3")`) recover most of it.
Estimating the SD of a heavy-tailed trait lowers coverage too, for gamma
and, with skewed traits, for beta even on a straight surface (90%), and
HC3 doesn't help there. Skew also biases beta a little. Refitting with a
transformed trait changes the scale selection is measured on, so report
it with the untransformed result. Collinearity is the largest variance
inflation factor from the linear gradient model. The rows-per-term check
counts the individuals, or for binary fitness the rarer outcome, per
term of the quadratic model, with ten as the working minimum. With more
than 2000 individuals Mardia's statistics use a random subsample of
2000, drawn the same way every run. Like the gradient models, the
logistic and count models fit an intercept per group, with unlabelled
rows as one more group. A group in which every individual survived, or
none did, is fitted at 0 or 1 by its intercept and is left out of the
separation check. For counts the note names the model behind each set of
p-values; the linear and quadratic models are checked for overdispersion
separately. If the `performance` package is installed its
heteroscedasticity and overdispersion tests are added next to the
package's own, along with an R squared for the gradient model.

## References

Lande, R. and Arnold, S. J. (1983) The measurement of selection on
correlated characters. Evolution 37, 1210-1226. Mardia, K. V. (1970)
Measures of multivariate skewness and kurtosis with applications.
Biometrika 57, 519-530. Morrissey, M. B. and Sakrejda, K. (2013)
Unification of regression-based methods for the analysis of natural
selection. Evolution 67, 2094-2100. Palacio, F. X., Ordano, M. and
Benitez-Vieyra, S. (2019) Measuring natural selection on multivariate
phenotypic traits: a protocol for verifiable and reproducible analyses
of natural selection. Israel Journal of Ecology and Evolution 65,
130-136.

## Examples

``` r
check_selection_assumptions(bumpus, "survival", c("total_length", "weight", "humerus"))
#> Assumption checks (binary fitness, n = 136)
#> 
#>  check                                    statistic p    
#>  Multivariate normality: Mardia skewness  1.05      0.006
#>  Multivariate normality: Mardia kurtosis    15      0.997
#>  Normality of total_length (Shapiro-Wilk) 0.978     0.029
#>  Normality of weight (Shapiro-Wilk)       0.97      0.004
#>  Normality of humerus (Shapiro-Wilk)      0.98      0.048
#>  Collinearity: largest VIF                1.75           
#>  Rarer outcome per quadratic term         7.11           
#>  Logistic model: separation               1.34           
#>  Model fit: R2_Tjur (performance)         0.24           
#>  note                                                      
#>  skewed                                                    
#>  no evidence against normality                             
#>  not normal                                                
#>  not normal                                                
#>  not normal                                                
#>  below 5                                                   
#>  64 of the rarer outcome for 9 terms; treat gamma with care
#>  no sign of separation                                     
#>                                                            
```
