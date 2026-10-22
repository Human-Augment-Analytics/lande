library(testthat)

# The app stores each heavy step under its data and arguments. These tests
# change one input at a time and check that a Run refits exactly the steps
# that take it, then that going back takes them all from the store. The
# plots are left undrawn: only the fits matter here.

app_steps <- c("report", "per_group", "assumptions", "canonical", "bootstrap", "band", "uni_land", "surface", "landscape")
plot_outputs <- c("hist_plot", "grad_plot", "uni_plot", "uni_land_plot", "surf_plot", "land_plot", "land3d",
                  "land3d_plotly", "group_plot")
app_inputs <- list(fitness = "survival", traits = c("weight", "total_length"), group = "sex", within_group = TRUE,
                   per_group = TRUE, ftype = "auto", surf_method = "gam", surf_grid = 20, grid_n = 15,
                   land_method = "lattice", sim_n = 100, n_boot = 20, seed = 1, surf_full = FALSE, surf_far = NA,
                   clamp = TRUE, canonical = TRUE, surf_k = NA, surf_bs = "tp", surf_sm = "REML", spline_bs = "cr",
                   spline_sm = "GCV.Cp", count_family = "poisson", uni_land = TRUE, spline_k = 5,
                   uni_trait = "weight", surf_x = "weight", surf_y = "total_length", classic = FALSE,
                   show_points = TRUE, show_groups = TRUE, point_alpha = 0.5, theme = "Viridis, orange and blue",
                   surf_unc = "none", show_opt = TRUE, show_mean = TRUE, show_support = FALSE, group_lines = TRUE,
                   land_theme = "Viridis, orange and blue", land_view = "2d")

# one session of the app: body() gets run_with(), which sets inputs, presses
# Run (after the k slider's 600 ms, unless `.wait = FALSE`), opens every
# result and presses the bootstrap button, and returns the steps fitted along
# the way. The band's refits go to two workers, which the app stops after the
# session if `last`.
app_session <- function(app, inputs, body, last = TRUE) {
  state <- new.env()
  state$runs <- 0
  old <- options(lande.workers = 2)
  on.exit(options(old), add = TRUE)
  shiny::testServer(app, {
    for (o in plot_outputs) session$defineOutput(o, NULL, o)
    # the test sees a copy of the server's variables, so the log is read
    # where the server keeps it
    log <- function() get("fits_made", envir = environment(stored))
    run_with <- function(..., .wait = TRUE) {
      n <- length(log())
      session$setInputs(...)
      if (.wait) session$elapse(700)
      state$runs <- state$runs + 1
      session$setInputs(run = state$runs)
      invisible(assump_checks())
      invisible(canon())
      session$setInputs(run_boot = state$runs)
      utils::tail(log(), length(log()) - n)
    }
    do.call(session$setInputs, inputs)
    state$value <- tryCatch(body(run_with),
                            finally = if (last) get("stop_band_cluster", envir = environment(stored))())
  })
  state$value
}

test_that("each input refits only the steps that take it, and going back refits nothing", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("fields")
  # half of the Bumpus birds, uploaded, to keep the fits quick, with made-up
  # counts of young as a second fitness
  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv), add = TRUE)
  half <- lande::bumpus[seq(1, nrow(lande::bumpus), 2), ]
  half$young <- (seq_len(nrow(half)) * 7) %% 5
  utils::write.csv(half, csv, row.names = FALSE)
  inputs <- c(list(dataset = "Upload CSV...", file = list(name = "half.csv", datapath = csv)), app_inputs)
  grad <- c("report", "per_group", "assumptions", "canonical", "bootstrap")
  band <- c("band", "uni_land")
  surf <- c("surface", "landscape")
  changes <- list(
    list(set = list(fitness = "head_length"), fits = app_steps),
    list(set = list(traits = c("weight", "total_length", "humerus")), fits = app_steps),
    list(set = list(group = "(none)"), fits = setdiff(app_steps, "per_group")),
    # without the group term each step but the surface (which still marks the
    # groups) and its landscape is the one just fitted with no group, so it
    # comes from the store
    list(set = list(within_group = FALSE), fits = surf),
    list(set = list(ftype = "continuous"), fits = c(grad, band)),
    list(set = list(seed = 2), fits = c("canonical", "bootstrap", band, "landscape")),
    list(set = list(n_boot = 30), fits = "bootstrap"),
    list(set = list(surf_method = "tps"), fits = surf),
    list(set = list(surf_grid = 30), fits = surf),
    list(set = list(surf_full = TRUE), fits = surf),
    list(set = list(surf_far = 0.2), fits = surf),
    list(set = list(surf_k = 8), fits = surf),
    list(set = list(surf_bs = "cr"), fits = surf),
    list(set = list(surf_sm = "ML"), fits = surf),
    list(set = list(clamp = FALSE), fits = c(surf, "uni_land")),
    # the count family reaches only count fits: not survival's, then the
    # band's once survival is fitted as a count (the surface still finds 0/1
    # and keeps its entry), then both on counts. Each pair's second change
    # finds the gradients stored by the first, which differs only in family.
    list(set = list(count_family = "nb"), fits = character()),
    list(set = list(ftype = "count", count_family = "quasipoisson"), fits = c(grad, band)),
    list(set = list(ftype = "count"), fits = band),
    list(set = list(fitness = "young", count_family = "quasipoisson"), fits = app_steps),
    list(set = list(fitness = "young"), fits = c(band, surf)),
    list(set = list(grid_n = 20), fits = c("uni_land", "landscape")),
    list(set = list(sim_n = 200), fits = c("uni_land", "landscape")),
    list(set = list(land_method = "convolve"), fits = c("uni_land", "landscape")),
    list(set = list(spline_bs = "tp"), fits = band),
    list(set = list(spline_sm = "REML"), fits = band),
    list(set = list(spline_k = 7), fits = band),
    list(set = list(uni_trait = "total_length"), fits = band),
    list(set = list(surf_x = "total_length", surf_y = "weight"), fits = surf),
    list(set = list(per_group = FALSE), fits = character()),
    list(set = list(canonical = FALSE), fits = character()),
    list(set = list(uni_land = FALSE), fits = character()),
    list(set = list(theme = "Magma, white and black", classic = TRUE, show_points = FALSE, point_alpha = 0.3,
                    surf_unc = "band", land_theme = "Greys, orange and blue", show_opt = FALSE, show_mean = FALSE,
                    show_groups = FALSE, group_lines = FALSE, show_support = TRUE, land_view = "3d"),
         fits = character())
  )
  app <- shiny::shinyAppDir(system.file("app", package = "lande"))
  app_session(app, inputs, function(run_with) {
    expect_setequal(run_with(), app_steps)
    for (ch in changes) {
      label <- paste(names(ch$set), collapse = ", ")
      got <- do.call(run_with, ch$set)
      expect(setequal(got, ch$fits), sprintf("changing %s refitted [%s], not [%s]", label,
                                             paste(got, collapse = ", "), paste(ch$fits, collapse = ", ")))
      back <- do.call(run_with, inputs[names(ch$set)])
      expect(length(back) == 0, sprintf("going back on %s refitted [%s]", label, paste(back, collapse = ", ")))
    }
    # a Run straight after the k slider moves fits the k it is on, and the
    # slider's own update 600 ms later refits nothing
    expect_setequal(run_with(spline_k = 9, .wait = FALSE), band)
    expect_length(run_with(), 0)
    expect_length(run_with(spline_k = inputs$spline_k, .wait = FALSE), 0)
  })
})

test_that("a bundled dataset's results are shared between sessions, an upload's are not", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  csv <- tempfile(fileext = ".csv")
  on.exit(unlink(csv), add = TRUE)
  utils::write.csv(lande::bumpus[seq(1, nrow(lande::bumpus), 2), ], csv, row.names = FALSE)
  quick <- utils::modifyList(app_inputs, list(per_group = FALSE, canonical = FALSE, uni_land = FALSE))
  bundled <- c(list(dataset = "Bumpus sparrows"), quick)
  uploaded <- c(list(dataset = "Upload CSV...", file = list(name = "half.csv", datapath = csv)), quick)
  fitted <- c("report", "assumptions", "bootstrap", "band", "surface", "landscape")
  app <- shiny::shinyAppDir(system.file("app", package = "lande"))
  run <- function(inputs, last = FALSE) app_session(app, inputs, function(run_with) run_with(), last = last)
  expect_setequal(run(bundled), fitted)
  expect_length(run(bundled), 0)
  expect_setequal(run(uploaded), fitted)
  expect_setequal(run(uploaded, last = TRUE), fitted)
})
