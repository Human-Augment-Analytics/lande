#' Bumpus house sparrow survival data
#'
#' Morphological measurements and overwinter survival of house sparrows
#' (\emph{Passer domesticus}) collected by Hermon Bumpus in 1898 after a severe
#' winter storm in Providence, Rhode Island.
#'
#' @format A data frame with 136 rows and 11 variables:
#' \describe{
#'   \item{sex}{Sex of the bird (\code{male} or \code{female}).}
#'   \item{survival}{Survival indicator: \code{1} if the bird survived, \code{0} otherwise.}
#'   \item{total_length}{Total length (mm).}
#'   \item{wingspread}{Alar extent / wingspread (mm).}
#'   \item{weight}{Body weight (g).}
#'   \item{head_length}{Length of head and beak (mm).}
#'   \item{humerus}{Length of humerus (in).}
#'   \item{femur}{Length of femur (in).}
#'   \item{tibiotarsus}{Length of tibiotarsus (in).}
#'   \item{skull_width}{Width of skull (in).}
#'   \item{sternum}{Length of sternum / keel (in).}
#' }
#' @source Bumpus, H. C. (1899). The elimination of the unfit as illustrated by
#'   the introduced sparrow, \emph{Passer domesticus}. \emph{Biological Lectures,
#'   Marine Biological Laboratory, Woods Hole}, 209-226.
"bumpus"

#' Pupfish hybrid survival in two San Salvador lakes
#'
#' Laboratory-reared F2 hybrids of the three pupfish species of San Salvador
#' Island, Bahamas, photographed, tagged and released into field enclosures in
#' Crescent Pond and Little Lake from March to June 2011, with
#' laboratory-reared fish of the three species measured the same way. Martin
#' (2016) analysed the high-density enclosures, \code{density == "H"}, as do
#' the examples and the app.
#'
#' The trait columns are Martin's scores for the 16 measurements of Martin and
#' Wainwright (2013), taken from the photographs. He log-transformed the
#' distances, regressed each trait on a size index (the first principal
#' component of seven size-related distances) and kept the residuals, fitting
#' the regression to the parental fish alone for \code{jaw}, \code{mouth},
#' \code{nasal} and \code{width}; \code{eye2}, \code{nose} and
#' \code{noseangle} were not size-corrected. Every trait is standardised to
#' mean 0 and standard deviation 1 within each lake, hybrids and parental fish
#' together.
#'
#' The column names are Martin's. The trait names below are the ones defined
#' in the supplement to Martin and Wainwright (2013), matched to the columns
#' through his size-correction script (in the Dryad archive) and their
#' discriminant loadings, with the name Martin (2016) uses in brackets where
#' it differs. The six functional traits of his Table 3 are \code{jaw},
#' \code{pmx}, \code{nose}, \code{noseangle}, \code{body} and \code{eye}. Their
#' loadings table swaps the labels of \code{snout} and \code{nasal}; the names
#' here follow the size-correction script.
#'
#' @format Data frames with 23 variables, 993 rows for Crescent Pond and 1062
#'   for Little Lake:
#' \describe{
#'   \item{lake}{\code{"CP"} (Crescent Pond) or \code{"LL"} (Little Lake).}
#'   \item{density}{The enclosure of an F2 hybrid, \code{"H"} (high density) or
#'     \code{"L"} (low density), or the species of a parental fish:
#'     \code{"norm"}, the generalist \emph{Cyprinodon variegatus};
#'     \code{"bozo"}, the molluscivore \emph{C. brontotheroides};
#'     \code{"bull"}, the scale-eater \emph{C. desquamator}.}
#'   \item{survival}{\code{1} if the hybrid survived the three months in the
#'     enclosure, \code{0} otherwise; \code{NA} for parental fish.}
#'   \item{ln.growth}{Log growth rate of survivors over the three months.
#'     Hybrids that died are coded \code{0}, so drop them before analysing
#'     growth; \code{NA} for parental fish.}
#'   \item{color}{Martin's plotting colour: orange for hybrids; blue, green
#'     and red for the generalist, molluscivore and scale-eater.}
#'   \item{d13C, d15N}{Carbon and nitrogen stable isotope ratios (per mil) of
#'     muscle from surviving hybrids; \code{NA} for the rest.}
#'   \item{jaw}{Lower jaw length, from the jaw joint to the tip of the
#'     dentary.}
#'   \item{eye}{Eye diameter, the mean of the major and minor axes of the iris
#'     (orbit diameter).}
#'   \item{eye2}{Eye roundness, the ratio of the major to the minor axis of the
#'     iris.}
#'   \item{pmx}{Craniofacial height, from the jaw joint to the tip of the
#'     premaxilla (upper jaw length).}
#'   \item{snout}{Lateral snout length, from the tip of the premaxilla to the
#'     front of the iris.}
#'   \item{body}{Dorsal to anal distance, between the first rays of the dorsal
#'     and anal fins (body depth).}
#'   \item{caudal}{Caudal peduncle height.}
#'   \item{nasal}{Dorsal snout length, from the front of the orbit to the tip
#'     of the maxilla in dorsal view.}
#'   \item{mouth}{Buccal width in dorsal view.}
#'   \item{width}{Head width across the opercula in dorsal view.}
#'   \item{boteyeangle}{Ventral orbit angle, at the tip of the premaxilla
#'     between the lower edge of the iris and the quadrate.}
#'   \item{topeyeangle}{Orbit angle, at the tip of the premaxilla between the
#'     upper and lower edges of the iris.}
#'   \item{nose}{Maxillary head protrusion (nasal protrusion).}
#'   \item{noseangle}{Maxillary head protrusion angle (nasal angle).}
#'   \item{adduct}{Preopercular height.}
#'   \item{SL}{Standard length.}
#' }
#' @source Martin, C. H. (2016) Context dependence in complex adaptive
#'   landscapes: frequency and trait-dependent selection surfaces within an
#'   adaptive radiation of Caribbean pupfishes. \emph{Evolution} 70,
#'   1265-1282. \doi{10.1111/evo.12932}. Data: \doi{10.5061/dryad.n3mj3}.
#' @references Martin, C. H. and Wainwright, P. C. (2013) Multiple fitness
#'   peaks on the adaptive landscape drive adaptive radiation in the wild.
#'   \emph{Science} 339, 208-211. \doi{10.1126/science.1227710}.
"crescent_pond_pupfish"

#' @rdname crescent_pond_pupfish
"little_lake_pupfish"

#' Medium ground finch survival from year to year
#'
#' Medium ground finches (\emph{Geospiza fortis}) marked and recaptured at El
#' Garrapatero, Santa Cruz Island, Galapagos. Each row is a bird in a year it
#' was seen, from 2004 to 2010, with whether it was seen again in any later
#' year (apparent survival). Built by \code{data-raw/finch_yearly.R}.
#'
#' @format A data frame with 1087 rows and 7 variables:
#' \describe{
#'   \item{band}{Band code of the bird.}
#'   \item{year}{Year the bird was seen.}
#'   \item{survived}{\code{1} if the bird was seen again in any later year
#'     up to 2018, the last in the file, \code{0} otherwise. Death cannot be
#'     told from emigration.}
#'   \item{beak_pc1}{Beak size: the first principal component of the three
#'     beak measurements over all the birds, signed so that larger beaks score
#'     higher.}
#'   \item{beak_length}{Median beak length (mm).}
#'   \item{beak_width}{Median beak width (mm).}
#'   \item{beak_depth}{Median beak depth (mm).}
#' }
#' @source Beausoleil, M.-O. et al. (2019) Temporally varying disruptive
#'   selection in the medium ground finch (\emph{Geospiza fortis}).
#'   \emph{Proceedings of the Royal Society B} 286, 20192290.
#'   \doi{10.1098/rspb.2019.2290}. Data: \doi{10.5061/dryad.zcrjdfn6q}.
"finch_yearly"

#' Darwin's finch community at El Garrapatero
#'
#' Four ground finch species at El Garrapatero, Santa Cruz Island, Galapagos,
#' with \emph{Geospiza fortis} split into its small and large beak morphs. One
#' row per bird, with its mean beak measurements and its apparent lifespan, the
#' fitness measure of Beausoleil et al. (2023). Birds first caught late in the
#' study had fewer years in which to be seen again. Built from the authors'
#' data deposit on Borealis (\doi{10.5683/SP3/0YIWSE}), released under CC BY
#' 4.0.
#'
#' @format A data frame with 3428 rows and 7 variables:
#' \describe{
#'   \item{band}{Band code of the bird.}
#'   \item{species}{A factor: \code{"fortis small"} and \code{"fortis large"}
#'     (\emph{G. fortis}), \code{"fuliginosa"} (\emph{G. fuliginosa}),
#'     \code{"magnirostris"} (\emph{G. magnirostris}) or \code{"scandens"}
#'     (\emph{G. scandens}).}
#'   \item{beak_length}{Mean beak length (mm).}
#'   \item{beak_depth}{Mean beak depth (mm).}
#'   \item{beak_width}{Mean beak width (mm).}
#'   \item{lifespan}{Apparent lifespan in years: the last year the bird was
#'     seen minus the first.}
#'   \item{first_year}{The year the bird was first caught.}
#' }
#' @source Beausoleil, M.-O. et al. (2023) The fitness landscape of a
#'   community of Darwin's finches. \emph{Evolution} 77, 2533-2546.
#'   \doi{10.1093/evolut/qpad160}. Data: \doi{10.5683/SP3/0YIWSE}.
"finch_community"
