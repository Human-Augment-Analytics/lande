# lande GUI

Shiny front-end for the package. Pick a dataset or upload a CSV (one row per
individual), choose the fitness column, the traits and optionally a group
column, then Run.

Tabs:

- **Data**: n, survival rate or mean fitness, group sizes, summary table,
  histograms, and the settings used (including the file name).
- **Selection gradients**: S, β ± SE and γ ± SE with p-values, correlational
  terms, a one-line summary per trait, a gradient plot, optional bootstrap
  intervals, CSV and PNG downloads, and the assumption checks (normality of
  the traits, VIF, rows per term, residuals or dispersion of the model).
- **Fitness functions**: cubic-spline fitness function per trait with a
  bootstrap band, and optionally the adaptive landscape for that trait.
- **Fitness surface**: two-trait surface with the individuals overlaid at an
  opacity you set, drawn only where there are data, in 2D or as a rotatable
  3D view. Six colour themes, each
  a surface palette with point colours that stand out against it. With a
  group, each group's mean and the highest point of the surface within its
  own range are marked, joined by a dashed line (switch in Advanced
  settings). Untick "Standardise within group" to treat the groups as one
  population, for several species on one surface.
- **Adaptive landscape**: 2D contour and rotatable 3D view, with the optimum
  and the current population mean, and its own colour theme for both views.
- **Groups**: gradients per group, table and plot.

Bootstrap, spline band and landscape simulation all run from the seed in
Advanced settings, so repeated runs match. Advanced settings also hold the
fitness type, the landscape method (lattice by default, which makes the same
draws as the original simulation; Gaussian convolution; or the simulation
itself), the spline and surface basis (cubic regression, thin plate or
P-spline) and smoothing criterion (GCV, REML or ML), the surface basis size
k, a switch to draw the surface over the whole grid rather than only where
there are data, a distance rule that blanks cells far from any individual,
and the switch for the line from each group mean to its peak.

Run shows the data and gradients first and fits the rest under the progress
bar. Each fit is kept, by its data and settings, so a Run after changing one
setting refits only what that setting feeds, and going back to earlier
settings refits nothing. The bundled datasets' fits are shared by every
session of the R process; an upload's stay in its session. Where the machine
has cores and memory to spare, the fitness function's 200 bootstrap refits
run on up to four worker processes (each about 230 MB), started at the first
Run and stopped with the app. They are not started on Windows, where the
workers' socket would bring up the firewall's prompt, or on shinyapps.io or
Posit Connect, unless `options(lande.workers = n)` with n of two or more asks
for n of them; `options(lande.workers = 0)` keeps the refits in the R session
everywhere. The results are the same either way. The band follows the k
slider once it has rested for 600 ms, and Run takes the value it is on.

## Run locally

```r
lande::run_app()          # from the installed package
shiny::runApp("inst/app")         # from a checkout, after devtools::load_all()
```

Bundled datasets: Bumpus sparrows, Crescent Pond and Little Lake pupfish (the
high-density enclosures Martin analysed), and the finch community of Beausoleil
et al. (2023), five groups with apparent lifespan as fitness, which opens with one
surface for all groups, quasi-Poisson counts and the 0.15 distance rule set. Each is a package
dataset with its own help page. Uploaded CSVs need numeric fitness and trait
columns.

The Data tab ends with the R calls that repeat the analysis with the current
settings, to read or download as a script.

## Optional packages

- `fields`: thin-plate spline surfaces, static 3D landscape and surface
- `plotly`: rotatable 3D landscape and surface
- `car`: Type III ANOVA tables

## GitHub Pages via Shinylive

```r
shinylive::export("app", "docs")
# commit docs/ and set Pages to deploy from /docs
```

Everything runs in the visitor's browser (first load fetches WebR and the
packages, roughly 20 to 40 MB). Keep the landscape grid and simulation sliders
modest. lande has no compiled code, so under WebR you load the CRAN
dependencies and source `R/*.R` directly.
