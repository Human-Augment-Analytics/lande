# Darwin's finch community at El Garrapatero

Four ground finch species at El Garrapatero, Santa Cruz Island,
Galapagos, with *Geospiza fortis* split into its small and large beak
morphs. One row per bird, with its mean beak measurements and its
apparent lifespan, the fitness measure of Beausoleil et al. (2023).
Birds first caught late in the study had fewer years in which to be seen
again. Built from the authors' data deposit on Borealis
([doi:10.5683/SP3/0YIWSE](https://doi.org/10.5683/SP3/0YIWSE) ),
released under CC BY 4.0.

## Usage

``` r
finch_community
```

## Format

A data frame with 3428 rows and 7 variables:

- band:

  Band code of the bird.

- species:

  A factor: `"fortis small"` and `"fortis large"` (*G. fortis*),
  `"fuliginosa"` (*G. fuliginosa*), `"magnirostris"` (*G. magnirostris*)
  or `"scandens"` (*G. scandens*).

- beak_length:

  Mean beak length (mm).

- beak_depth:

  Mean beak depth (mm).

- beak_width:

  Mean beak width (mm).

- lifespan:

  Apparent lifespan in years: the last year the bird was seen minus the
  first.

- first_year:

  The year the bird was first caught.

## Source

Beausoleil, M.-O. et al. (2023) The fitness landscape of a community of
Darwin's finches. *Evolution* 77, 2533-2546.
[doi:10.1093/evolut/qpad160](https://doi.org/10.1093/evolut/qpad160) .
Data: [doi:10.5683/SP3/0YIWSE](https://doi.org/10.5683/SP3/0YIWSE) .
