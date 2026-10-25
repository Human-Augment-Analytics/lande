# Medium ground finch survival from year to year

Medium ground finches (*Geospiza fortis*) marked and recaptured at El
Garrapatero, Santa Cruz Island, Galapagos. Each row is a bird in a year
it was seen, from 2004 to 2010, with whether it was seen again in any
later year (apparent survival). Built by `data-raw/finch_yearly.R` in
the GitHub repository.

## Usage

``` r
finch_yearly
```

## Format

A data frame with 1087 rows and 7 variables:

- band:

  Band code of the bird.

- year:

  Year the bird was seen.

- survived:

  `1` if the bird was seen again in any later year up to 2018, the last
  in the file, `0` otherwise. Death cannot be told from emigration.

- beak_pc1:

  Beak size: the first principal component of the three beak
  measurements over all the birds, signed so that larger beaks score
  higher.

- beak_length:

  Median beak length (mm).

- beak_width:

  Median beak width (mm).

- beak_depth:

  Median beak depth (mm).

## Source

Beausoleil, M.-O. et al. (2019) Temporally varying disruptive selection
in the medium ground finch (*Geospiza fortis*). *Proceedings of the
Royal Society B* 286, 20192290.
[doi:10.1098/rspb.2019.2290](https://doi.org/10.1098/rspb.2019.2290) .
Data:
[doi:10.5061/dryad.zcrjdfn6q](https://doi.org/10.5061/dryad.zcrjdfn6q) .
