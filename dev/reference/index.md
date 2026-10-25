# Package index

## Prepare

Standardise the traits, relativise fitness and check the data

- [`prepare_selection_data()`](https://human-augment-analytics.github.io/lande/dev/reference/prepare_selection_data.md)
  : Prepare data for selection analysis
- [`detect_family()`](https://human-augment-analytics.github.io/lande/dev/reference/detect_family.md)
  : Automatically detect fitness type and GLM family
- [`check_selection_assumptions()`](https://human-augment-analytics.github.io/lande/dev/reference/check_selection_assumptions.md)
  : Check the assumptions behind a selection analysis

## Differentials and gradients

Selection differentials, linear and quadratic gradients, their intervals
and canonical axes

- [`selection_differential()`](https://human-augment-analytics.github.io/lande/dev/reference/selection_differential.md)
  : Calculate selection differential (S)
- [`selection_coefficients()`](https://human-augment-analytics.github.io/lande/dev/reference/selection_coefficients.md)
  : Calculate selection coefficients
- [`analyze_linear_selection()`](https://human-augment-analytics.github.io/lande/dev/reference/analyze_linear_selection.md)
  : Analyze linear selection gradients (beta)
- [`analyze_nonlinear_selection()`](https://human-augment-analytics.github.io/lande/dev/reference/analyze_nonlinear_selection.md)
  : Analyze nonlinear selection gradients (gamma)
- [`analyze_disruptive_selection()`](https://human-augment-analytics.github.io/lande/dev/reference/analyze_disruptive_selection.md)
  : Disruptive or stabilising selection on one trait
- [`extract_linear_coefficients()`](https://human-augment-analytics.github.io/lande/dev/reference/extract_linear_coefficients.md)
  : Extract linear selection coefficients
- [`extract_quadratic_coefficients()`](https://human-augment-analytics.github.io/lande/dev/reference/extract_quadratic_coefficients.md)
  : Extract quadratic selection coefficients
- [`extract_interaction_coefficients()`](https://human-augment-analytics.github.io/lande/dev/reference/extract_interaction_coefficients.md)
  : Extract correlational selection coefficients
- [`bootstrap_selection()`](https://human-augment-analytics.github.io/lande/dev/reference/bootstrap_selection.md)
  : Bootstrap selection gradients
- [`canonical_analysis()`](https://human-augment-analytics.github.io/lande/dev/reference/canonical_analysis.md)
  : Canonical analysis of the gamma matrix
- [`plot_canonical_axes()`](https://human-augment-analytics.github.io/lande/dev/reference/plot_canonical_axes.md)
  : Fitness along the canonical axes
- [`selection_report()`](https://human-augment-analytics.github.io/lande/dev/reference/selection_report.md)
  : Standardised selection analysis report

## Fitness functions and surfaces

Spline fitness functions for one trait and fitness surfaces for two

- [`univariate_spline()`](https://human-augment-analytics.github.io/lande/dev/reference/univariate_spline.md)
  : Estimate univariate correlated fitness function
- [`plot_univariate_fitness()`](https://human-augment-analytics.github.io/lande/dev/reference/plot_univariate_fitness.md)
  : Plot Univariate Correlated Fitness Function
- [`correlated_fitness_surface()`](https://human-augment-analytics.github.io/lande/dev/reference/correlated_fitness_surface.md)
  : Calculate the Correlated Fitness Surface
- [`plot_correlated_fitness()`](https://human-augment-analytics.github.io/lande/dev/reference/plot_correlated_fitness.md)
  : Plot Correlated Fitness Surface
- [`plot_correlated_fitness_enhanced()`](https://human-augment-analytics.github.io/lande/dev/reference/plot_correlated_fitness_enhanced.md)
  : Plot Enhanced Correlated Fitness Surface
- [`plot_correlated_fitness_3d()`](https://human-augment-analytics.github.io/lande/dev/reference/plot_correlated_fitness_3d.md)
  : Plot Correlated Fitness Surface (3D Perspective)
- [`peak_difference()`](https://human-augment-analytics.github.io/lande/dev/reference/peak_difference.md)
  : Compare the fitted fitness at two points of a surface

## Adaptive landscape

Mean fitness against the population mean, and how it compares with the
individual surface

- [`adaptive_landscape()`](https://human-augment-analytics.github.io/lande/dev/reference/adaptive_landscape.md)
  : Calculate Adaptive Landscape
- [`plot_adaptive_landscape()`](https://human-augment-analytics.github.io/lande/dev/reference/plot_adaptive_landscape.md)
  : Plot Adaptive Landscape
- [`plot_adaptive_landscape_3d()`](https://human-augment-analytics.github.io/lande/dev/reference/plot_adaptive_landscape_3d.md)
  : Plot Adaptive Landscape (3D Perspective)
- [`compare_fitness_surfaces_data()`](https://human-augment-analytics.github.io/lande/dev/reference/compare_fitness_surfaces_data.md)
  : Compare Correlated Fitness Surface vs Adaptive Landscape
- [`plot_fitness_surfaces_comparison()`](https://human-augment-analytics.github.io/lande/dev/reference/plot_fitness_surfaces_comparison.md)
  : Plot Fitness Surfaces Comparison

## Over time

One fitness function or surface per period

- [`temporal_landscape()`](https://human-augment-analytics.github.io/lande/dev/reference/temporal_landscape.md)
  : Fitness functions and adaptive landscapes over time
- [`plot_temporal_landscape()`](https://human-augment-analytics.github.io/lande/dev/reference/plot_temporal_landscape.md)
  : Plot fitness functions and landscapes over time

## App

The whole workflow in a browser

- [`run_app()`](https://human-augment-analytics.github.io/lande/dev/reference/run_app.md)
  : Open the Shiny app
- [`lande`](https://human-augment-analytics.github.io/lande/dev/reference/lande-package.md)
  [`lande-package`](https://human-augment-analytics.github.io/lande/dev/reference/lande-package.md)
  : lande: Standardised Analyses of Natural Selection Data

## Datasets

Example datasets included with lande

- [`bumpus`](https://human-augment-analytics.github.io/lande/dev/reference/bumpus.md)
  : Bumpus house sparrow survival data
- [`crescent_pond_pupfish`](https://human-augment-analytics.github.io/lande/dev/reference/crescent_pond_pupfish.md)
  [`little_lake_pupfish`](https://human-augment-analytics.github.io/lande/dev/reference/crescent_pond_pupfish.md)
  : Pupfish hybrid survival in two San Salvador lakes
- [`finch_yearly`](https://human-augment-analytics.github.io/lande/dev/reference/finch_yearly.md)
  : Medium ground finch survival from year to year
- [`finch_community`](https://human-augment-analytics.github.io/lande/dev/reference/finch_community.md)
  : Darwin's finch community at El Garrapatero

## Print methods

- [`print(`*`<selection_report>`*`)`](https://human-augment-analytics.github.io/lande/dev/reference/print.selection_report.md)
  : Print a selection report
- [`print(`*`<adaptive_landscape>`*`)`](https://human-augment-analytics.github.io/lande/dev/reference/print.adaptive_landscape.md)
  : Print method for adaptive landscape
- [`print(`*`<selection_assumptions>`*`)`](https://human-augment-analytics.github.io/lande/dev/reference/print.selection_assumptions.md)
  : Print the assumption checks
