# Plot Correlated Fitness Surface (3D Perspective)

Draws the fitted fitness of a two-trait surface in perspective. Cells
blanked by `mask` or `too_far` are left out, so only the part of the
surface with data under it is drawn. Each grid cell is outlined, so a
surface fitted with `grid_n = 80` shows as 79 x 79 small squares.

## Usage

``` r
plot_correlated_fitness_3d(
  tps,
  trait_cols,
  theta = -30,
  phi = 30,
  color_palette = NULL,
  show_peaks = TRUE,
  border = "grey30",
  lwd = 0.5,
  shade = 0.15,
  ...
)
```

## Arguments

- tps:

  Output list from
  [`correlated_fitness_surface()`](https://human-augment-analytics.github.io/lande/dev/reference/correlated_fitness_surface.md).

- trait_cols:

  Character vector of length 2 specifying the trait column names.

- theta:

  Numeric azimuthal viewing angle. Default is -30.

- phi:

  Numeric colatitude viewing angle. Default is 30.

- color_palette:

  Colours for the surface: a vector of colours, or `"jet"` for the
  rainbow ramp (dark blue through cyan, yellow and red to dark red) used
  in many published fitness surfaces. Defaults to viridis. Jet is harder
  to read for colour-blind viewers and its bright bands can look like
  steps that are not in the surface, so viridis is the safer choice for
  most figures.

- show_peaks:

  Logical; mark the interior peaks of the surface (the rows of
  `tps$peaks` with `interior` `TRUE`) with gold diamonds. Maxima at the
  edge of the data are not marked. Default is `TRUE`.

- border:

  Colour of the lines outlining each grid cell. Default is `"grey30"`, a
  dark grey that keeps the squares clear without darkening the colours
  as black lines do; `NA` draws the surface without them. With more than
  100 grid values along a side the default lines are faded, more as the
  grid gets finer, so they do not hide the colours; give `border` to
  keep them at full strength.

- lwd:

  Width of the cell outlines, which also sets the width of the box
  lines. Default is 0.5.

- shade:

  Strength of the shading that shows the relief, from 0 (flat colours,
  no shading) towards 1 (darker slopes), passed to
  [`graphics::persp()`](https://rdrr.io/r/graphics/persp.html). Default
  is 0.15, light enough that the colours stay close to the colour bar
  while the hills still stand out. Use `NA` to switch the shading off.

- ...:

  Additional arguments passed to
  [`fields::drape.plot()`](https://rdrr.io/pkg/fields/man/drape.plot.html).

## Value

The viewing transformation matrix from
[`fields::drape.plot()`](https://rdrr.io/pkg/fields/man/drape.plot.html),
invisibly, which
[`grDevices::trans3d()`](https://rdrr.io/r/grDevices/trans3d.html) uses
to add points.

## Examples

``` r
prep <- prepare_selection_data(bumpus, "survival", c("total_length", "weight"))
surf <- correlated_fitness_surface(prep, "survival", c("total_length", "weight"), grid_n = 30)
#> Data type: binary; method: gam; n = 136; k = 29
#> GAM fitting with 136 observations
#>   Trying formula: main
#> Success with formula: main
#> Predictions range: 0.0472 to 0.703
#> Masked 424 of 900 grid points outside the data
plot_correlated_fitness_3d(surf, c("total_length", "weight"))
```
