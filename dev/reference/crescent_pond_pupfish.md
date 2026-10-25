# Pupfish hybrid survival in two San Salvador lakes

Laboratory-reared F2 hybrids of the three pupfish species of San
Salvador Island, Bahamas, photographed, tagged and released into field
enclosures in Crescent Pond and Little Lake from March to June 2011,
with laboratory-reared fish of the three species measured the same way.
Martin (2016) analysed the high-density enclosures, `density == "H"`, as
do the examples and the app.

## Usage

``` r
crescent_pond_pupfish

little_lake_pupfish
```

## Format

Data frames with 23 variables, 993 rows for Crescent Pond and 1062 for
Little Lake:

- lake:

  `"CP"` (Crescent Pond) or `"LL"` (Little Lake).

- density:

  The enclosure of an F2 hybrid, `"H"` (high density) or `"L"` (low
  density), or the species of a parental fish: `"norm"`, the generalist
  *Cyprinodon variegatus*; `"bozo"`, the molluscivore *C.
  brontotheroides*; `"bull"`, the scale-eater *C. desquamator*.

- survival:

  `1` if the hybrid survived the three months in the enclosure, `0`
  otherwise; `NA` for parental fish.

- ln.growth:

  Log growth rate of survivors over the three months. Hybrids that died
  are coded `0`, so drop them before analysing growth; `NA` for parental
  fish.

- color:

  Martin's plotting colour: orange for hybrids; blue, green and red for
  the generalist, molluscivore and scale-eater.

- d13C, d15N:

  Carbon and nitrogen stable isotope ratios (per mil) of muscle from
  surviving hybrids; `NA` for the rest.

- jaw:

  Lower jaw length, from the jaw joint to the tip of the dentary.

- eye:

  Eye diameter, the mean of the major and minor axes of the iris (orbit
  diameter).

- eye2:

  Eye roundness, the ratio of the major to the minor axis of the iris.

- pmx:

  Craniofacial height, from the jaw joint to the tip of the premaxilla
  (upper jaw length).

- snout:

  Lateral snout length, from the tip of the premaxilla to the front of
  the iris.

- body:

  Dorsal to anal distance, between the first rays of the dorsal and anal
  fins (body depth).

- caudal:

  Caudal peduncle height.

- nasal:

  Dorsal snout length, from the front of the orbit to the tip of the
  maxilla in dorsal view.

- mouth:

  Buccal width in dorsal view.

- width:

  Head width across the opercula in dorsal view.

- boteyeangle:

  Ventral orbit angle, at the tip of the premaxilla between the lower
  edge of the iris and the quadrate.

- topeyeangle:

  Orbit angle, at the tip of the premaxilla between the upper and lower
  edges of the iris.

- nose:

  Maxillary head protrusion (nasal protrusion).

- noseangle:

  Maxillary head protrusion angle (nasal angle).

- adduct:

  Preopercular height.

- SL:

  Standard length.

An object of class `data.frame` with 1062 rows and 23 columns.

## Source

Martin, C. H. (2016) Context dependence in complex adaptive landscapes:
frequency and trait-dependent selection surfaces within an adaptive
radiation of Caribbean pupfishes. *Evolution* 70, 1265-1282.
[doi:10.1111/evo.12932](https://doi.org/10.1111/evo.12932) . Data:
[doi:10.5061/dryad.n3mj3](https://doi.org/10.5061/dryad.n3mj3) .

## Details

The trait columns are Martin's scores for the 16 measurements of Martin
and Wainwright (2013), taken from the photographs. He log-transformed
the distances, regressed each trait on a size index (the first principal
component of seven size-related distances) and kept the residuals,
fitting the regression to the parental fish alone for `jaw`, `mouth`,
`nasal` and `width`; `eye2`, `nose` and `noseangle` were not
size-corrected. Every trait is standardised to mean 0 and standard
deviation 1 within each lake, hybrids and parental fish together.

The column names are Martin's. The trait names below are the ones
defined in the supplement to Martin and Wainwright (2013), matched to
the columns through his size-correction script (in the Dryad archive)
and their discriminant loadings, with the name Martin (2016) uses in
brackets where it differs. The six functional traits of his Table 3 are
`jaw`, `pmx`, `nose`, `noseangle`, `body` and `eye`. Their loadings
table swaps the labels of `snout` and `nasal`; the names here follow the
size-correction script.

## References

Martin, C. H. and Wainwright, P. C. (2013) Multiple fitness peaks on the
adaptive landscape drive adaptive radiation in the wild. *Science* 339,
208-211.
[doi:10.1126/science.1227710](https://doi.org/10.1126/science.1227710) .
