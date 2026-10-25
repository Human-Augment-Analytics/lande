# Plot Univariate Correlated Fitness Function

Visualizes the univariate fitness function estimated via spline.

## Usage

``` r
plot_univariate_fitness(
  uni,
  trait_col,
  title = NULL,
  show_points = TRUE,
  point_alpha = 0.25,
  ribbon_fill = NA,
  ribbon_linetype = "dashed",
  classic_plot = FALSE,
  ...
)
```

## Arguments

- uni:

  Output object of class `"univariate_fitness"` from
  [`univariate_spline()`](https://human-augment-analytics.github.io/lande/dev/reference/univariate_spline.md).

- trait_col:

  Character string specifying the trait column name.

- title:

  Optional character string for the plot title.

- show_points:

  Logical indicating whether to show the original data points. Default
  is `TRUE`.

- point_alpha:

  Numeric value for point transparency. Default is 0.25.

- ribbon_fill:

  Color for the confidence interval ribbon. Default is `NA`
  (transparent).

- ribbon_linetype:

  Linetype for the confidence interval bounds. Default is `"dashed"`.

- classic_plot:

  Logical; if `TRUE`, draw the plain published style (a single dashed
  confidence outline, solid fit line, raw points, and a minimal theme).
  Default is `FALSE`.

- ...:

  Additional arguments passed to
  [`ggplot2::labs()`](https://ggplot2.tidyverse.org/reference/labs.html).

## Value

A `ggplot` object representing the univariate fitness function.

## Examples

``` r
prep <- prepare_selection_data(bumpus, "survival", "total_length")
uni <- univariate_spline(prep, "survival", "total_length")
#> Fitness type detected: binary
plot_univariate_fitness(uni, "total_length")

plot_univariate_fitness(uni, "total_length", classic_plot = TRUE)
```
