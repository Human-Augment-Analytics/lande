# Builds inst/extdata/finch_community.csv and data/finch_community.rda: the
# finches of El Garrapatero, four species with Geospiza fortis split into small
# and large beak morphs, from the public data of Beausoleil et al. (2023), one
# row per bird with its mean beak measurements and its apparent lifespan, the
# last year it was seen minus the first (their fitness measure, mxcpois in their
# file). Source file kept in validation/data:
# it is data/bird.data.RData from the authors' Borealis deposit,
# https://doi.org/10.5683/SP3/0YIWSE, released under CC BY 4.0 (the same file
# is in their GPL-3 code repository). Cite the paper and the dataset.
# Run from the package root.
load("validation/data/beausoleil2023_bird.data.RData")
d <- bird.data[complete.cases(bird.data[, c("avg.mbl", "avg.mbd", "avg.mbw", "mxcpois", "sp2")]), ]
years <- 2003:2020
seen <- as.matrix(d[, paste0("y.", years)]) == 1
out <- data.frame(
  band = d$BANDFINAL,
  species = d$sp2,
  beak_length = d$avg.mbl,
  beak_depth = d$avg.mbd,
  beak_width = d$avg.mbw,
  lifespan = as.integer(d$mxcpois),
  # birds first caught late had fewer years in which to be seen again
  first_year = apply(seen, 1, function(r) min(years[which(r)]))
)
rownames(out) <- NULL
write.csv(out, "inst/extdata/finch_community.csv", row.names = FALSE)
finch_community <- out
save(finch_community, file = "data/finch_community.rda", compress = "xz")
