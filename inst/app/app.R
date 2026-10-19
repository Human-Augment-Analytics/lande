# ---------------------------------------------------------------------------
# lande: interactive selection analysis
#
# Run locally:   shiny::runApp("app")
# Static export (GitHub Pages): shinylive::export("app", "docs"); serve docs/
#
# The core analysis needs only mgcv and MASS. The thin-plate spline surface,
# the 3D landscape (fields) and Type III ANOVA (car) degrade to a message
# when those packages are missing.
# ---------------------------------------------------------------------------

library(shiny)
library(lande)
library(ggplot2)

`%||%` <- function(a, b) if (is.null(a)) b else a
# spelled out so rsconnect picks them up
has_fields <- requireNamespace("fields", quietly = TRUE)
has_plotly <- requireNamespace("plotly", quietly = TRUE)
has_car <- requireNamespace("car", quietly = TRUE)
has_performance <- requireNamespace("performance", quietly = TRUE)

# fixed seed for the random steps; restores the RNG afterwards
with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv())
          else suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE)
  set.seed(seed)
  expr
}

# run a fit and keep its warnings, to show beside the results. The VIF warning
# is left out: it comes from the quadratic model, where squares and products
# are collinear by construction, and the assumption table has the linear VIF.
with_warnings <- function(expr) {
  warned <- character()
  value <- withCallingHandlers(suppressMessages(expr), warning = function(w) {
    warned <<- c(warned, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = grep("VIF above 5", warned, value = TRUE, invert = TRUE))
}

# download names start with the dataset: a short name for the bundled ones,
# the file name for an upload
STEMS <- c("Bumpus sparrows" = "bumpus", "Crescent Pond pupfish" = "crescent_pond_pupfish",
           "Little Lake pupfish" = "little_lake_pupfish", "Finch community (five groups)" = "finch_community")
file_stem <- function(source) {
  bundled <- sub(" \\(bundled\\)$", "", source)
  if (bundled %in% names(STEMS)) return(STEMS[[bundled]])
  stem <- gsub("^_+|_+$", "", tolower(gsub("[^A-Za-z0-9]+", "_", sub("\\.[^.]*$", "", basename(source)))))
  if (nzchar(stem)) stem else "data"
}

# Martin analysed the high-density enclosures only
high_density <- function(d) d[d$density == "H", ]

load_dataset <- function(name) {
  switch(name,
    "Bumpus sparrows" = list(
      data = lande::bumpus, fitness = "survival",
      traits = c("weight", "total_length"), group = "sex"),
    "Crescent Pond pupfish" = list(
      data = high_density(lande::crescent_pond_pupfish),
      fitness = "survival", traits = c("jaw", "body"), group = NULL),
    "Little Lake pupfish" = list(
      data = high_density(lande::little_lake_pupfish),
      fitness = "survival", traits = c("jaw", "body"), group = NULL),
    # five groups on one surface: standardised together, blank far from any bird
    "Finch community (five groups)" = list(
      data = lande::finch_community,
      fitness = "lifespan", traits = c("beak_length", "beak_depth"), group = "species",
      within_group = FALSE, too_far = 0.15, count_family = "quasipoisson")
  )
}

DATASETS <- c("Bumpus sparrows", "Crescent Pond pupfish", "Little Lake pupfish",
              "Finch community (five groups)", "Upload CSV...")
numeric_cols <- function(df) names(df)[vapply(df, is.numeric, logical(1))]

# --- formatting helpers -----------------------------------------------------
stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
}
fmt <- function(x, d = 3) ifelse(is.na(x), "", formatC(x, digits = d, format = "f"))
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "< 0.001", paste("=", formatC(p, digits = 3, format = "f"))))
pm <- function(est, se) paste0(fmt(est), " ± ", fmt(se))
sq <- function(t) paste0(t, "²")

# colour themes for the surface and landscape plots: a palette for the surface
# and two point colours that stand out against it (perished/survived, or the
# low and high ends of a continuous fitness scale)
THEMES <- list(
  "Viridis, orange and blue" = list(fill = function(...) ggplot2::scale_fill_viridis_d(option = "D", ...),
                                    ramp = function(n) viridisLite::viridis(n, option = "D"),
                                    points = c("#D55E00", "#0072B2")),
  "Magma, white and black"   = list(fill = function(...) ggplot2::scale_fill_viridis_d(option = "A", ...),
                                    ramp = function(n) viridisLite::viridis(n, option = "A"),
                                    points = c("white", "black")),
  "Plasma, teal and white"   = list(fill = function(...) ggplot2::scale_fill_viridis_d(option = "C", ...),
                                    ramp = function(n) viridisLite::viridis(n, option = "C"),
                                    points = c("#004D40", "white")),
  "Cividis, red and blue"    = list(fill = function(...) ggplot2::scale_fill_viridis_d(option = "E", ...),
                                    ramp = function(n) viridisLite::viridis(n, option = "E"),
                                    points = c("#B2182B", "#4393C3")),
  "Greys, orange and blue"   = list(fill = function(...) ggplot2::scale_fill_grey(start = 0.92, end = 0.25, ...),
                                    ramp = function(n) grDevices::grey.colors(n, start = 0.92, end = 0.25),
                                    points = c("#D55E00", "#0072B2"))
)
theme_of <- function(theme) if (!is.null(theme) && theme %in% names(THEMES)) THEMES[[theme]] else THEMES[[1]]
# shorter legend labels for the binned surface and landscape fills, "0.064–0.13"
# for "(0.0642, 0.1283]", with just enough digits to keep the bins apart
short_bins <- function(lv) {
  txt <- regmatches(lv, gregexpr("-?[0-9]*\\.?[0-9]+(e[-+]?[0-9]+)?", lv))
  if (!length(lv) || !all(lengths(txt) == 2)) return(lv)
  lo <- as.numeric(vapply(txt, `[`, "", 1)); hi_txt <- vapply(txt, `[`, "", 2); hi <- as.numeric(hi_txt)
  # ggplot prints the two ends of a bin separately, so one boundary can come as
  # 0.147 and 0.1467: count it once where the two agree at the printed digits
  for (i in seq_len(length(lv) - 1)) {
    if (!grepl("e", hi_txt[i])) {
      dec <- nchar(sub("^[^.]*\\.?", "", hi_txt[i]))
      if (round(lo[i + 1], dec) == hi[i]) hi[i] <- lo[i + 1]
    }
  }
  b <- sort(unique(c(lo, hi)))
  for (d in 2:6) {
    f <- trimws(formatC(b, digits = d, format = "fg"))
    if (!anyDuplicated(f)) break
  }
  f <- sub("^-", "−", f)
  sep <- if (any(b < 0)) " to " else "–"
  paste0(f[match(lo, b)], sep, f[match(hi, b)])
}
apply_theme <- function(p, theme, binary = NULL, fill_name = NULL, point_name = "Fitness") {
  th <- theme_of(theme)
  p <- p + th$fill(name = fill_name, labels = short_bins)
  if (isTRUE(binary)) {
    p <- p + ggplot2::scale_colour_manual(values = c("0" = th$points[1], "1" = th$points[2]),
                                          labels = c("0" = "Perished", "1" = "Survived"), name = "Outcome")
  } else if (isFALSE(binary)) {
    p <- p + ggplot2::scale_colour_gradient(low = th$points[1], high = th$points[2], name = point_name)
  }
  p
}

# one row per trait: S, beta, gamma
gradient_table <- function(r, traits) {
  pick <- function(type, term) {
    i <- which(r$Type == type & r$Term == term)
    if (length(i)) r[i[1], ] else data.frame(Estimate = NA, Std_Error = NA, P_Value = NA)
  }
  rows <- lapply(traits, function(t) {
    S <- pick("Differential", t); L <- pick("Linear", t); Q <- pick("Quadratic", sq(t))
    data.frame(
      Trait = t,
      S = fmt(S$Estimate),
      beta = pm(L$Estimate, L$Std_Error),
      p_beta = paste(sub("^= ", "", fmt_p(L$P_Value)), stars(L$P_Value)),
      gamma = pm(Q$Estimate, Q$Std_Error),
      p_gamma = paste(sub("^= ", "", fmt_p(Q$P_Value)), stars(Q$P_Value)),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  names(out) <- c("Trait", "S", "β ± SE", "p (linear term)", "γ ± SE", "p (squared term)")
  out
}

correlational_table <- function(r) {
  cc <- r[r$Type == "Correlational", , drop = FALSE]
  if (!nrow(cc)) return(NULL)
  data.frame(
    Traits = cc$Term,
    `γij ± SE` = pm(cc$Estimate, cc$Std_Error),
    p = paste(sub("^= ", "", fmt_p(cc$P_Value)), stars(cc$P_Value)),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

# the R calls that repeat the analysis outside the app, with the current settings
r_value <- function(x) {
  if (is.null(x)) "NULL" else if (is.character(x)) {
    if (length(x) == 1) sprintf('"%s"', x) else sprintf("c(%s)", paste(sprintf('"%s"', x), collapse = ", "))
  } else if (is.logical(x)) as.character(x) else format(x)
}
# one call as text, wrapped at 80 characters under its bracket
r_call <- function(fn, ..., assign = NULL) {
  args <- Filter(Negate(is.null), list(...))
  nm <- names(args)
  if (is.null(nm)) nm <- rep("", length(args))
  vals <- vapply(args, function(a) as.character(a)[1], "")
  parts <- ifelse(nzchar(nm), paste(nm, "=", vals), vals)
  head <- paste0(if (!is.null(assign)) paste0(assign, " <- "), fn, "(")
  one <- paste0(head, paste(parts, collapse = ", "), ")")
  if (nchar(one) <= 80) return(one)
  pad <- strrep(" ", nchar(head))
  out <- head
  cur <- nchar(head)
  for (i in seq_along(parts)) {
    piece <- paste0(parts[i], if (i < length(parts)) "," else ")")
    if (cur + nchar(piece) + 1 > 80 && cur > nchar(head)) {
      out <- paste0(out, "\n", pad)
      cur <- nchar(pad)
    } else if (i > 1) {
      out <- paste0(out, " ")
      cur <- cur + 1
    }
    out <- paste0(out, piece)
    cur <- cur + nchar(piece)
  }
  out
}
DATA_CODE <- list(
  "Bumpus sparrows" = "dat <- bumpus",
  "Crescent Pond pupfish" = c("dat <- crescent_pond_pupfish",
                              'dat <- dat[dat$density == "H", ]  # the enclosures Martin analysed'),
  "Little Lake pupfish" = c("dat <- little_lake_pupfish",
                            'dat <- dat[dat$density == "H", ]  # the enclosures Martin analysed'),
  "Finch community (five groups)" = "dat <- finch_community"
)
r_code <- function(s, dataset, file_name, uni_trait, spline_k, surf_traits, n_boot, uncertainty, canonical,
                   uni_land = FALSE) {
  fit <- r_value(s$fit); grp <- r_value(s$group_model); type <- r_value(s$ftype)
  cf <- if (identical(s$ftype, "count")) r_value(s$count_family)
  # the surface finds counts in the data even when the fitness type is set otherwise,
  # so it takes the app's count family in any case
  cf_surf <- if (!identical(s$count_family, "poisson")) r_value(s$count_family)
  load <- if (dataset %in% names(DATA_CODE)) DATA_CODE[[dataset]] else
    sprintf('dat <- read.csv("%s", na.strings = c("NA", ""))', file_name %||% "your_file.csv")
  lines <- c(
    "library(lande)", "", load,
    paste("traits <-", r_value(s$traits)),
    # the app keeps only rows with fitness and every trait; unlabelled rows are one more group
    sprintf("dat <- dat[complete.cases(dat[, c(%s, traits)]), ]", fit),
    "", "# differentials and gradients, with the checks to report beside them",
    r_call("selection_report", "dat", fit, "traits", fitness_type = type, group = grp),
    if (s$per_group) r_call("selection_coefficients", "dat", fit, "traits", fitness_type = type, group = r_value(s$group), return_grouped = "TRUE"),
    r_call("check_selection_assumptions", "dat", fit, "traits", fitness_type = type, group = grp),
    sprintf("set.seed(%d)", s$seed),
    r_call("bootstrap_selection", "dat", fit, "traits", fitness_type = type, group = grp, n_boot = n_boot),
    if (canonical && length(s$traits) > 1) c(sprintf("set.seed(%d)", s$seed),
                                             r_call("canonical_analysis", "dat", fit, "traits", fitness_type = type, group = grp)),
    "", "# fitness function",
    r_call("prepare_selection_data", "dat", fit, "traits", group = grp, na_action = '"none"', assign = "prep"),
    sprintf("set.seed(%d)", s$seed),
    r_call("univariate_spline", "prep", fit, r_value(uni_trait), fitness_type = type, group = grp, k = spline_k,
           bs = r_value(s$spline_bs), smoothing = r_value(s$spline_sm), bootstrap = "TRUE", n_boot = 200,
           count_family = cf, assign = "uni"),
    r_call("plot_univariate_fitness", "uni", r_value(uni_trait)),
    if (uni_land) c(sprintf("set.seed(%d)", s$seed),
                    r_call("adaptive_landscape", "prep", "uni$model", r_value(uni_trait), group_col = grp, grid_n = s$grid_n,
                           simulation_n = s$sim_n, clamp = r_value(s$clamp), assign = "land1"),
                    r_call("plot_adaptive_landscape", "land1", r_value(uni_trait)))
  )
  if (length(surf_traits) == 2 && surf_traits[1] != surf_traits[2]) {
    tr <- r_value(surf_traits)
    lines <- c(lines, "", "# fitness surface and adaptive landscape",
      r_call("correlated_fitness_surface", "prep", fit, tr, method = r_value(s$surf_method), grid_n = s$surf_grid,
             mask = r_value(!s$surf_full), too_far = r_value(s$surf_far), group = r_value(s$group),
             group_effect = r_value(s$within_group), k = r_value(s$surf_k), bs = r_value(s$surf_bs),
             smoothing = r_value(s$surf_sm), clamp = r_value(s$clamp), count_family = cf_surf, assign = "surf"),
      r_call("plot_correlated_fitness", "surf", tr, uncertainty = r_value(uncertainty)),
      sprintf("set.seed(%d)", s$seed),
      r_call("adaptive_landscape", "prep", "surf$model", tr, group_col = grp, grid_n = s$grid_n,
             simulation_n = s$sim_n, clamp = r_value(s$clamp), assign = "land"),
      r_call("plot_adaptive_landscape", "land", tr))
  }
  lines
}

# short reading of the gradient table
interpret <- function(r, traits, ftype, n, group) {
  get <- function(type, term) {
    i <- which(r$Type == type & r$Term == term)
    if (length(i)) r[i[1], ] else NULL
  }
  sig <- function(x) !is.null(x) && !is.na(x$P_Value) && x$P_Value < 0.05
  est <- function(x, sym) sprintf("%s = %s, p %s", sym, pm(x$Estimate, x$Std_Error), fmt_p(x$P_Value))
  lines <- vapply(traits, function(t) {
    L <- get("Linear", t); Q <- get("Quadratic", sq(t))
    if (is.null(L)) return(sprintf("%s: not estimated.", t))
    if (!sig(L) && !sig(Q)) {
      return(sprintf("%s: no detectable selection (%s%s).", t, est(L, "β"),
                     if (is.null(Q)) "" else paste0("; ", est(Q, "γ"))))
    }
    d <- if (sig(L)) sprintf("selection for %s values (%s)", if (L$Estimate > 0) "larger" else "smaller", est(L, "β"))
         else sprintf("no clear evidence of directional selection (%s)", est(L, "β"))
    q <- if (is.null(Q)) NULL
         else if (sig(Q)) sprintf("%s curvature, consistent with %s selection (%s)",
                                  if (Q$Estimate < 0) "negative" else "positive",
                                  if (Q$Estimate < 0) "stabilising" else "disruptive", est(Q, "γ"))
         else sprintf("no clear evidence of curvature (%s)", est(Q, "γ"))
    sprintf("%s: %s.", t, paste(c(d, q), collapse = "; "))
  }, character(1))
  cc <- r[r$Type == "Correlational", , drop = FALSE]
  if (nrow(cc)) {
    cl <- vapply(seq_len(nrow(cc)), function(i) {
      x <- cc[i, ]
      if (!is.na(x$P_Value) && x$P_Value < 0.05) {
        sprintf("%s: correlational selection for %s combinations (%s).", x$Term,
                if (x$Estimate > 0) "matching" else "opposing", est(x, "γij"))
      } else sprintf("%s: no clear evidence of correlational selection (%s).", x$Term, est(x, "γij"))
    }, character(1))
    lines <- c(lines, cl)
  }
  within <- if (is.null(group)) "" else sprintf(" within %s", group)
  # name the count model that was fitted; the quadratic one can differ
  plain <- function(f) c("poisson(log)" = "Poisson", "negative binomial" = "negative binomial")[f]
  pm <- attr(r, "p_model")
  count_model <- if (length(pm) == 2 && all(pm %in% c("poisson(log)", "negative binomial"))) {
    if (pm[1] == pm[2]) sprintf("a %s model", plain(pm[1])) else sprintf("a %s model for β and a %s model for γ and γij", plain(pm[1]), plain(pm[2]))
  } else "a Poisson or negative binomial model"
  basis <- switch(ftype,
    binary = " Binary fitness: p-values are for the terms of a logistic model, not β and γ themselves.",
    count = sprintf(" Count fitness: p-values are for the terms of %s, not β and γ themselves.", count_model),
    ""
  )
  c(lines, "", sprintf("n = %d. Traits in SD units, fitness relative to the mean%s.%s γ is curvature; check the fitness function for a peak.", n, within, basis))
}

# traits are standardised before every fit, so the plot axes are in SD units
sd_lab <- function(t, mean = FALSE) paste0(if (mean) "mean " else "", t, " (SD)")
# links open in a new tab, so the session and its results stay open
doc_link <- function(text, href) tags$a(href = href, target = "_blank", rel = "noopener", text)
# a result's heading with its downloads on the right
res_head <- function(title, ...) div(class = "res-head", h4(class = "sec", title), div(class = "res-dl", ...))
dl_btn <- function(id, label, what) downloadButton(id, label, title = what, `aria-label` = what)

# --- UI ---------------------------------------------------------------------
ui <- fluidPage(
  title = "lande",
  tags$head(tags$style(HTML(
    "body{background:#F7F8F5;color:#24332D}
     a,a:focus{color:#356653}
     a:hover{color:#2B5444}
     .well{background:#FFFFFF;border-color:#DEE5DF;box-shadow:none}
     .nav-tabs{border-bottom-color:#DEE5DF}
     .nav-tabs>li>a{color:#356653}
     .nav-tabs>li.active>a,.nav-tabs>li.active>a:hover,.nav-tabs>li.active>a:focus{color:#24332D;border-color:#DEE5DF #DEE5DF transparent}
     .tab-content{background:#FFFFFF;border:1px solid #DEE5DF;border-top:none;border-radius:0 0 4px 4px;padding:0 16px 16px}
     .btn-primary,.btn-primary:focus{background:#356653;border-color:#2B5444}
     .btn-primary:hover,.btn-primary:active,.btn-primary.active{background:#2B5444;border-color:#22433A}
     .form-control,.selectize-input{border-color:#7D8F85}
     .form-control:focus,.selectize-input.focus{border-color:#356653;box-shadow:0 0 0 2px rgba(53,102,83,.25)}
     .selectize-dropdown .active{background:#E3ECE6;color:#24332D}
     .selectize-control.multi .selectize-input>div{background:#E3ECE6;color:#24332D}
     .irs--shiny .irs-bar{background:#356653;border-color:#356653}
     .irs--shiny .irs-single,.irs--shiny .irs-from,.irs--shiny .irs-to{background:#356653}
     input[type=checkbox],input[type=radio]{accent-color:#356653}
     .progress-bar{background-color:#356653}
     .help-note{color:#4F5D57;font-size:13px;margin-top:6px}
     .help-note.warn{color:#8a4b00}
     .interp{background:#EEF3EF;border-left:3px solid #356653;padding:10px 14px;margin:10px 0;font-size:14px}
     .interp .shiny-text-output{white-space:pre-wrap;padding:0;margin:0}
     .interp ul{margin:0;padding-left:20px}
     .interp li{margin:5px 0}
     .interp p.foot{margin:10px 0 0;color:#4F5D57;font-size:13px}
     details summary{cursor:pointer;font-weight:600;margin:10px 0}
     h4.sec{margin-top:18px}
     .app-header{display:flex;align-items:center;flex-wrap:wrap;gap:6px 14px;padding:10px 20px 12px}
     .app-header img{height:52px;width:auto}
     .app-header h1{font-size:24px;font-weight:600;margin:0}
     .app-header .app-sub{color:#4F5D57;font-size:14px}
     .app-header .app-links{margin-left:auto;display:flex;flex-wrap:wrap;gap:16px;font-size:14px}
     .intro{background:#EEF3EF;border-left:3px solid #356653;padding:10px 14px;margin:0 0 12px;font-size:14px}
     .res-head{display:flex;align-items:flex-end;justify-content:space-between;flex-wrap:wrap;gap:8px;margin:18px 0 8px}
     .res-head h4.sec{margin:0}
     .res-head .res-dl{display:flex;align-items:center;flex-wrap:wrap;gap:6px}
     .res-head .res-dl .form-group{margin:0}
     .res-head .res-dl .radio-inline{padding-top:0}
     @media (max-width:767px){#land_plot{height:400px !important}}"
  ))),
  tags$head(tags$link(rel = "icon", type = "image/png", href = "favicon.png")),
  # drop the connection after ten minutes without a click or keypress, so a tab
  # left open does not keep the server instance awake
  tags$head(tags$script(HTML("
    (function () {
      var idleMinutes = 10;
      var timer;
      function reset() {
        clearTimeout(timer);
        timer = setTimeout(function () {
          if (window.Shiny && Shiny.shinyapp && Shiny.shinyapp.$socket) Shiny.shinyapp.$socket.close();
          document.body.innerHTML = '<p style=\"font-family: sans-serif; margin: 40px\">Disconnected after ' + idleMinutes + ' minutes without activity. Reload the page to start again.</p>';
        }, idleMinutes * 60 * 1000);
      }
      ['mousemove', 'keydown', 'click', 'touchstart', 'scroll'].forEach(function (e) {
        document.addEventListener(e, reset, true);
      });
      reset();
    })();
  "))),
  div(class = "app-header",
      tags$img(src = "logo.png", alt = ""),
      div(tags$h1("lande"), div(class = "app-sub", "Evolutionary selection analysis")),
      div(class = "app-links",
          doc_link("Guide", "https://human-augment-analytics.github.io/lande/articles/evolutionary-selection-analysis.html"),
          doc_link("Reference", "https://human-augment-analytics.github.io/lande/reference/"),
          doc_link("Report a problem", "https://github.com/Human-Augment-Analytics/lande/issues"))),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Data"),
      selectInput("dataset", "Dataset", DATASETS),
      conditionalPanel("input.dataset == 'Upload CSV...'",
        fileInput("file", "CSV file", accept = ".csv")),
      h4("Design"),
      selectInput("fitness", "Fitness", NULL),
      div(class = "help-note", textOutput("fitness_hint")),
      selectizeInput("traits", "Traits", NULL, multiple = TRUE),
      selectInput("group", "Group (optional)", "(none)"),
      checkboxInput("within_group", "Standardise within group, group as a fixed effect", TRUE),
      div(class = "help-note",
          "Untick to treat all groups as one population; the group then only marks means and peaks on the surface, as for several species on one surface."),
      checkboxInput("per_group", "Estimate per group as well", FALSE),
      actionButton("run", "Run analysis", class = "btn-primary", width = "100%"),
      tags$details(
        tags$summary("Advanced settings"),
        selectInput("ftype", "Fitness type", c("detect automatically" = "auto", "binary (0/1)" = "binary", "count" = "count", "continuous" = "continuous")),
        radioButtons("surf_method", "Fitness-surface method",
                     c("GAM (mgcv)" = "gam", "Thin-plate spline (fields)" = "tps"), inline = TRUE),
        sliderInput("surf_grid", "Surface grid", 20, 80, 50, step = 10),
        sliderInput("grid_n", "Landscape grid", 15, 60, 40, step = 5),
        sliderInput("sim_n", "Simulated individuals per grid point", 100, 1000, 300, step = 100),
        numericInput("n_boot", "Bootstrap resamples", 500, 100, 5000, 100),
        numericInput("seed", "Random seed", 1, 1, 1e6, 1),
        checkboxInput("surf_full", "Draw the surface beyond the data", FALSE),
        numericInput("surf_far", "Blank surface cells farther than this share of the axis range from any individual (blank: off)", NA, 0.02, 1, 0.01),
        checkboxInput("group_lines", "Join each group mean to its peak on the surface", TRUE),
        checkboxInput("clamp", "Keep thin-plate fitness within the range of the fitness type", TRUE),
        checkboxInput("canonical", "Canonical analysis of γ on the gradients tab", FALSE),
        numericInput("surf_k", "Surface basis size k (blank: from the data)", NA, 5, 60, 1),
        selectInput("surf_bs", "Surface basis (GAM)", c("thin plate" = "tp", "cubic regression" = "cr", "P-spline" = "ps")),
        selectInput("surf_sm", "Surface smoothing (GAM)", c("REML" = "REML", "GCV / UBRE" = "GCV.Cp", "ML" = "ML")),
        selectInput("spline_bs", "Fitness function basis", c("cubic regression" = "cr", "thin plate" = "tp", "P-spline" = "ps")),
        selectInput("spline_sm", "Fitness function smoothing", c("GCV / UBRE" = "GCV.Cp", "REML" = "REML", "ML" = "ML")),
        selectInput("count_family", "Count fitness family (GAM)", c("Poisson" = "poisson", "quasi-Poisson" = "quasipoisson", "negative binomial" = "nb"))
      )
    ),
    mainPanel(
      width = 9,
      conditionalPanel("input.run == 0", div(class = "intro", textOutput("intro", inline = TRUE))),
      tabsetPanel(id = "tabs",
        tabPanel("Data",
          br(), verbatimTextOutput("data_summary"),
          h4(class = "sec", "Traits and fitness"), tableOutput("trait_summary"),
          plotOutput("hist_plot", height = "300px"),
          h4(class = "sec", "Settings"), verbatimTextOutput("settings"),
          res_head("R code", dl_btn("dl_code", "R script", "Download the R script")),
          div(class = "help-note", "The calls that repeat this analysis in R with the settings above."),
          verbatimTextOutput("r_code")),
        tabPanel("Selection gradients",
          res_head("Selection differentials and gradients",
                   dl_btn("dl_report", "CSV", "Download the gradients as CSV")),
          tableOutput("grad_table"),
          uiOutput("corr_table_ui"),
          uiOutput("fit_warnings"),
          uiOutput("canon_ui"),
          div(class = "help-note", "S total selection; β directional; γ quadratic curvature (negative is consistent with stabilising selection, positive with disruptive); γij correlational. * p < 0.05, ** < 0.01, *** < 0.001."),
          h4(class = "sec", "Summary"),
          uiOutput("interpretation"),
          res_head("Gradient plot", dl_btn("dl_gradplot", "PNG", "Download the gradient plot as PNG")),
          plotOutput("grad_plot", height = "360px"),
          actionButton("run_boot", "Add bootstrap confidence intervals"),
          div(class = "help-note", textOutput("boot_note")),
          res_head("Assumption checks", dl_btn("dl_assump", "CSV", "Download the assumption checks as CSV")),
          tableOutput("assump_table"),
          div(class = "help-note", "Normality of the traits (Mardia, Shapiro-Wilk), the largest VIF, individuals per quadratic term, and the residuals or dispersion of the gradient model. For reporting only. With the performance package installed its tests are added.")),
        tabPanel("Fitness functions",
          br(),
          fluidRow(
            column(4, selectInput("uni_trait", "Trait", NULL)),
            column(4, sliderInput("spline_k", "Spline flexibility (k)", 3, 12, 5, step = 1)),
            column(4, br(), checkboxInput("classic", "Classic style", FALSE))),
          res_head("Fitness function", dl_btn("dl_uniplot", "PNG", "Download the fitness function as PNG")),
          plotOutput("uni_plot", height = "440px"),
          div(class = "help-note", "Cubic spline (Schluter 1988) with bootstrap 95% band."),
          br(),
          checkboxInput("uni_land", "Adaptive landscape for this trait", FALSE),
          conditionalPanel("input.uni_land",
            res_head("Adaptive landscape for this trait",
                     dl_btn("dl_unilandplot", "PNG", "Download the one-trait landscape as PNG")),
            plotOutput("uni_land_plot", height = "400px"),
            div(class = "help-note", "Mean fitness of a population against its mean for this trait, with the fitness function dashed."))),
        tabPanel("Fitness surface",
          br(),
          fluidRow(
            column(3, selectInput("surf_x", "Trait on x", NULL)),
            column(3, selectInput("surf_y", "Trait on y", NULL)),
            column(3, br(), checkboxInput("show_points", "Show individuals", TRUE)),
            column(3, br(), checkboxInput("show_groups", "Group means and peaks", TRUE))),
          fluidRow(
            column(5, conditionalPanel("input.show_points",
              sliderInput("point_alpha", "Opacity of the individuals", 0.05, 1, 0.5, step = 0.05))),
            column(4, selectInput("theme", "Colours", names(THEMES))),
            column(3, selectInput("surf_unc", "Uncertainty (GAM)",
                                  c("none" = "none", "standard error lines" = "se", "lower, fit and upper" = "band")))),
          res_head("Fitness surface", dl_btn("dl_surfplot", "PNG", "Download the fitness surface as PNG")),
          plotOutput("surf_plot", height = "500px"),
          div(class = "help-note", "A group's peak is the highest point of the surface within that group's own range: filled when it is a peak of the surface, open when the surface keeps rising past the group's range or the edge of the data.")),
        tabPanel("Adaptive landscape",
          br(), div(class = "interp", textOutput("opt_txt")),
          fluidRow(
            column(2, br(), checkboxInput("show_opt", "Optimum", TRUE)),
            column(3, br(), checkboxInput("show_mean", "Current mean", TRUE)),
            column(3, conditionalPanel("input.land_view == '2d'", br(),
              checkboxInput("show_support", "Shade where over half the simulation is outside the data", FALSE))),
            column(4, selectInput("land_theme", "Colours", names(THEMES)))),
          res_head("Adaptive landscape",
                   radioButtons("land_view", tags$span(class = "sr-only", "Landscape view"),
                                c("2D" = "2d", "3D" = "3d"), inline = TRUE),
                   conditionalPanel("input.land_view == '2d'",
                     dl_btn("dl_landplot", "PNG", "Download the 2D landscape as PNG"))),
          conditionalPanel("input.land_view == '2d'",
            plotOutput("land_plot", height = "560px"),
            div(class = "help-note", "Mean fitness for each possible population mean. Diamond: optimum; circle: current mean.")),
          conditionalPanel("input.land_view == '3d'",
            uiOutput("land3d_ui"),
            div(class = "help-note", "Mean fitness for each possible population mean. Drag to rotate. Diamond: optimum; circle: current mean."))),
        tabPanel("Groups",
          br(), uiOutput("groups_ui"))
      )
    )
  )
)

# --- server -----------------------------------------------------------------
server <- function(input, output, session) {

  current <- reactive({
    if (input$dataset == "Upload CSV...") {
      req(input$file)
      # a blank cell is missing, as a spreadsheet writes it, not a group named ""
      c(list(data = utils::read.csv(input$file$datapath, na.strings = c("NA", "")), fitness = NULL, traits = NULL, group = NULL),
        source = input$file$name)
    } else c(load_dataset(input$dataset), source = paste(input$dataset, "(bundled)"))
  })

  # Defaults for an uploaded file: a survival-like or 0/1 column as fitness, a
  # year or site column as group, and measured traits (not ids or years) as traits
  observeEvent(current(), {
    d <- current(); nm <- names(d$data); nums <- numeric_cols(d$data)
    looks <- function(cols, pattern) cols[grepl(pattern, cols, ignore.case = TRUE)]
    n_unique <- function(v) length(unique(d$data[[v]][!is.na(d$data[[v]])]))
    is_01 <- function(v) all(d$data[[v]][!is.na(d$data[[v]])] %in% c(0, 1))
    id_like <- looks(nm, "^id$|band|ring|tag|^year$|site|sex|pop|species|lake|group|cohort|season|plot")
    fit_default <- d$fitness %||% c(
      looks(nums, "surv|fitness|fecund|offspring|seed|growth|^w$"),
      nums[vapply(nums, is_01, logical(1))],
      setdiff(nums, id_like), nums
    )[1]
    grp_default <- d$group %||% {
      g <- looks(nm, "^year$|site|sex|pop|species|lake|group|cohort|season")
      g <- g[vapply(g, function(v) n_unique(v) %in% 2:30, logical(1))]
      if (length(g)) g[1] else "(none)"
    }
    t_choices <- setdiff(nums, fit_default)
    t_default <- d$traits %||% {
      t <- setdiff(t_choices, c(id_like, grp_default))
      t <- t[vapply(t, function(v) n_unique(v) > 5, logical(1))]
      if (!length(t)) t <- t_choices
      t[seq_len(min(2, length(t)))]
    }
    updateSelectInput(session, "fitness", choices = nm, selected = fit_default)
    updateSelectizeInput(session, "traits", choices = t_choices, selected = t_default)
    updateSelectInput(session, "group", choices = c("(none)", nm), selected = grp_default)
    # dataset presets for the group handling and the distance rule
    updateCheckboxInput(session, "within_group", value = !isFALSE(d$within_group))
    updateNumericInput(session, "surf_far", value = d$too_far %||% NA)
    updateSelectInput(session, "count_family", selected = d$count_family %||% "poisson")
  })

  output$fitness_hint <- renderText({
    d <- current()$data; f <- input$fitness
    if (is.null(f) || !f %in% names(d) || !is.numeric(d[[f]])) return("")
    fam <- detect_family(d[[f]])
    if (fam$type == "binary") sprintf("Survival (0/1): %.0f%% of %d survived.", 100 * mean(d[[f]], na.rm = TRUE), sum(!is.na(d[[f]])))
    else sprintf("%s: mean %.2f, SD %.2f.", tools::toTitleCase(fam$type), mean(d[[f]], na.rm = TRUE), stats::sd(d[[f]], na.rm = TRUE))
  })

  # everything shared across tabs, computed once at Run
  setup <- eventReactive(input$run, {
    d <- current()$data; fit <- input$fitness; traits <- input$traits
    grp <- if (input$group == "(none)") NULL else input$group
    validate(
      need(length(traits) >= 1, "Pick at least one trait."),
      need(!fit %in% traits, "Fitness column can't also be a trait."),
      need(all(c(fit, traits) %in% names(d)), "Columns not found."),
      need(is.numeric(d[[fit]]) || is.logical(d[[fit]]), "Fitness column must be numeric (0/1 survival or a numeric measure)."),
      need(all(vapply(d[traits], is.numeric, logical(1))), "Trait columns must be numeric."),
      need(all(vapply(d[traits], function(x) length(unique(x[!is.na(x)])) > 2, logical(1))),
           "Each trait needs more than two distinct values; a 0/1 column cannot be a trait.")
    )
    # rows need fitness and every trait; the package keeps unlabelled rows as one more group
    d <- d[stats::complete.cases(d[, c(fit, traits)]), ]
    validate(need(nrow(d) >= 20, "Need at least 20 individuals with complete data."))
    resolved <- if (input$ftype == "auto") detect_family(d[[fit]])$type else input$ftype
    if (!resolved %in% c("binary", "count")) resolved <- "continuous"
    # with the switch off the groups are one population: pooled standardisation,
    # no group term, the group only marks means and peaks on the surface
    within <- !is.null(grp) && !isFALSE(input$within_group)
    grp_model <- if (within) grp else NULL
    prep <- suppressWarnings(suppressMessages(
      prepare_selection_data(d, fit, traits, standardize = TRUE, group = grp_model, add_relative = TRUE, na_action = "none")))
    fitted <- with_warnings(
      selection_report(d, fit, traits, fitness_type = resolved, standardize = TRUE, group = grp_model))
    report <- fitted$value
    warned <- fitted$warnings
    grouped <- NULL
    if (!is.null(grp) && isTRUE(input$per_group)) {
      per <- tryCatch(with_warnings(
        selection_coefficients(d, fit, traits, fitness_type = resolved, standardize = TRUE,
                               group = grp, return_grouped = TRUE)),
        error = function(e) list(value = NULL, warnings = paste("Per-group fits failed:", conditionMessage(e))))
      grouped <- per$value
      if (length(per$warnings)) warned <- c(warned, paste("Per group:", per$warnings))
    }
    keep_or <- function(x, default) if (isTRUE(x %in% traits)) x else default
    updateSelectInput(session, "uni_trait", choices = traits, selected = keep_or(input$uni_trait, traits[1]))
    updateSelectInput(session, "surf_x", choices = traits, selected = keep_or(input$surf_x, traits[1]))
    updateSelectInput(session, "surf_y", choices = traits, selected = keep_or(input$surf_y, traits[min(2, length(traits))]))
    list(d = d, prep = prep, fit = fit, traits = traits, group = grp, ftype = resolved,
         report = report, grouped = grouped, warnings = unique(warned), source = current()$source,
         # what was run, for the script, even if the inputs change afterwards
         dataset = input$dataset, file_name = if (!is.null(input$file)) input$file$name,
         ftype_how = if (input$ftype == "auto") "detected" else "set by hand",
         per_group = isTRUE(input$per_group),
         seed = if (is.numeric(input$seed) && !is.na(input$seed)) round(input$seed) else 1,
         surf_method = input$surf_method, surf_grid = input$surf_grid, surf_full = isTRUE(input$surf_full),
         surf_k = if (is.numeric(input$surf_k) && !is.na(input$surf_k)) round(input$surf_k) else NULL,
         surf_far = if (is.numeric(input$surf_far) && !is.na(input$surf_far) && input$surf_far > 0) input$surf_far else NULL,
         clamp = !isFALSE(input$clamp),
         within_group = within, group_model = grp_model,
         surf_bs = input$surf_bs %||% "tp", surf_sm = input$surf_sm %||% "REML",
         spline_bs = input$spline_bs %||% "cr", spline_sm = input$spline_sm %||% "GCV.Cp",
         count_family = input$count_family %||% "poisson",
         grid_n = input$grid_n, sim_n = input$sim_n)
  })

  boot_val <- reactiveVal(NULL)
  # Do the heavy work once at Run, so the tabs only draw afterwards
  observeEvent(input$run, {
    withProgress(message = "Running the analysis", value = 0, {
      incProgress(0.05, detail = "gradients")
      s <- try(setup(), silent = TRUE)
      req(!inherits(s, "try-error"))
      uni_choice(if (isTRUE(input$uni_trait %in% s$traits)) input$uni_trait else s$traits[1])
      surf_choice(pick_surface(s, input$surf_x, input$surf_y))
      incProgress(0.35, detail = "fitness function")
      try(uni_fit(), silent = TRUE)
      if (isTRUE(input$uni_land)) try(uni_land(), silent = TRUE)
      incProgress(0.35, detail = "fitness surface and landscape")
      if (length(s$traits) >= 2) try(surfaces(), silent = TRUE)
      incProgress(0.25, detail = "done")
    })
  }, priority = 100)
  observeEvent(input$run, boot_val(NULL))
  observeEvent(input$run_boot, {
    s <- setup()
    validate(need(is.numeric(input$n_boot) && !is.na(input$n_boot) && input$n_boot >= 2,
                  "Resamples must be a whole number of at least 2."))
    withProgress(message = "Bootstrapping...", {
      b <- with_seed(s$seed, suppressWarnings(suppressMessages(
        bootstrap_selection(s$d, s$fit, s$traits, fitness_type = s$ftype, standardize = TRUE,
                            group = s$group_model, n_boot = round(input$n_boot)))))
    })
    # the number drawn, for the script; attr(b, "n_boot") counts the usable ones
    attr(b, "n_drawn") <- round(input$n_boot)
    boot_val(b)
  })

  # ---- Before the first run ----
  output$intro <- renderText({
    what <- "Estimates selection differentials and gradients, fitness functions, fitness surfaces and adaptive landscapes."
    if (identical(input$dataset, "Upload CSV...")) {
      paste(what, "Upload a CSV with one row per individual, pick the fitness and traits, then press Run analysis.")
    } else {
      paste(what, sprintf("%s is loaded: press Run analysis, or pick another dataset or upload a CSV.", input$dataset))
    }
  })

  # ---- Data tab ----
  output$data_summary <- renderText({
    s <- setup(); f <- detect_family(s$d[[s$fit]])
    fit_line <- if (s$ftype == "binary") {
      sprintf("Fitness: %s, binary. %d of %d survived (%.0f%%).", s$fit, sum(s$d[[s$fit]] == 1), nrow(s$d), 100 * mean(s$d[[s$fit]]))
    } else {
      sprintf("Fitness: %s, %s. Mean %.3f, SD %.3f.", s$fit, f$type, mean(s$d[[s$fit]]), stats::sd(s$d[[s$fit]]))
    }
    grp_line <- if (is.null(s$group)) "No grouping: one population." else {
      tab <- table(s$d[[s$group]])
      sprintf("Grouped by %s: %s.", s$group, paste(sprintf("%s (n = %d)", names(tab), as.integer(tab)), collapse = ", "))
    }
    paste(sprintf("%d individuals with complete data.", nrow(s$d)), fit_line, grp_line,
          sprintf("Traits: %s.", paste(s$traits, collapse = ", ")), sep = "\n")
  })
  output$trait_summary <- renderTable({
    s <- setup()
    vars <- c(s$fit, s$traits)
    data.frame(
      Variable = vars,
      Mean = vapply(vars, function(v) mean(s$d[[v]]), 1),
      SD = vapply(vars, function(v) stats::sd(s$d[[v]]), 1),
      Min = vapply(vars, function(v) min(s$d[[v]]), 1),
      Max = vapply(vars, function(v) max(s$d[[v]]), 1),
      check.names = FALSE)
  }, digits = 3)
  output$hist_plot <- renderPlot({
    s <- setup()
    long <- do.call(rbind, lapply(c(s$fit, s$traits), function(v) data.frame(variable = v, value = s$d[[v]])))
    long$variable <- factor(long$variable, levels = c(s$fit, s$traits))
    ggplot(long, aes(value)) + geom_histogram(bins = 25, fill = "#2c7fb8", colour = "white") +
      facet_wrap(~variable, scales = "free", nrow = 1) + theme_bw() + labs(x = NULL, y = "Individuals")
  })

  # ---- Selection gradients tab ----
  output$grad_table <- renderTable({ s <- setup(); gradient_table(s$report, s$traits) }, align = "lrrrrr")
  assump_checks <- reactive({
    s <- setup()
    suppressWarnings(check_selection_assumptions(s$d, s$fit, s$traits, fitness_type = s$ftype, group = s$group_model))
  })
  output$assump_table <- renderTable({
    chk <- assump_checks()
    data.frame(
      Check = chk$check,
      Statistic = ifelse(is.na(chk$statistic), "", formatC(chk$statistic, digits = 3, format = "g")),
      p = ifelse(is.na(chk$p_value), "", ifelse(chk$p_value < 0.001, "< 0.001", formatC(chk$p_value, digits = 3, format = "f"))),
      Note = chk$note, check.names = FALSE, stringsAsFactors = FALSE)
  }, align = "lrrl")
  output$dl_assump <- downloadHandler(
    filename = function() paste0(dl_name("assumption_checks"), ".csv"),
    content = function(f) {
      chk <- assump_checks()
      utils::write.csv(data.frame(check = chk$check, statistic = chk$statistic, p_value = chk$p_value,
                                  note = chk$note), f, row.names = FALSE)
    })
  output$corr_table_ui <- renderUI({
    s <- setup(); if (is.null(correlational_table(s$report))) return(NULL)
    tagList(h4(class = "sec", "Correlational selection"), tableOutput("corr_table"))
  })
  output$corr_table <- renderTable({ s <- setup(); correlational_table(s$report) }, align = "lrr")
  output$fit_warnings <- renderUI({
    w <- setup()$warnings
    if (!length(w)) return(NULL)
    div(class = "help-note warn", paste("Warnings from the fit:", paste(w, collapse = "; ")))
  })

  # canonical axes of gamma, only when asked for in Advanced settings
  canon <- reactive({
    s <- setup()
    if (!isTRUE(input$canonical) || length(s$traits) < 2) return(NULL)
    # the p-values come from shuffling fitness, so from the stored seed
    tryCatch(with_seed(s$seed, suppressWarnings(suppressMessages(
      canonical_analysis(s$d, s$fit, s$traits, fitness_type = s$ftype, standardize = TRUE, group = s$group_model)))),
      error = function(e) NULL)
  })
  output$canon_ui <- renderUI({
    if (is.null(canon())) return(NULL)
    tagList(h4(class = "sec", "Canonical axes of γ"), tableOutput("canon_table"),
            div(class = "help-note", "λ is the curvature along each axis (negative is consistent with stabilising selection, positive with disruptive), and θ the directional selection along it. The p-values come from 999 shuffles of fitness among individuals (Reynolds et al. 2010), which allow for the axes being estimated from these data. The largest curvatures are still overestimated, and the SE treats the axes as known."))
  })
  output$canon_table <- renderTable({
    ca <- canon(); if (is.null(ca)) return(NULL)
    a <- ca$axes
    out <- data.frame(Axis = a$axis, `λ ± SE` = pm(a$lambda, a$se),
                      p = paste(sub("^= ", "", fmt_p(a$p_value)), stars(a$p_value)), `θ` = fmt(a$theta),
                      check.names = FALSE, stringsAsFactors = FALSE)
    load <- as.data.frame(t(round(ca$M, 2)), check.names = FALSE)
    cbind(out, load)
  }, align = "l")
  output$interpretation <- renderUI({
    # pass the group only if it was in the model, so the summary says "within" only then
    s <- setup(); lines <- interpret(s$report, s$traits, s$ftype, nrow(s$d), s$group_model)
    sep <- which(lines == "")[1]
    bullets <- lines[seq_len(sep - 1)]; foot <- lines[(sep + 1):length(lines)]
    items <- lapply(bullets, function(l) {
      t <- s$traits[startsWith(l, paste0(s$traits, ": "))]
      if (length(t) == 1) tags$li(HTML(paste0("<b>", htmltools::htmlEscape(t), "</b>", htmltools::htmlEscape(substring(l, nchar(t) + 1))))) else tags$li(l)
    })
    div(class = "interp", tags$ul(items), tags$p(class = "foot", paste(foot, collapse = " ")))
  })
  grad_plot_obj <- reactive({
    s <- setup(); r <- s$report[s$report$Type != "Differential", ]
    b <- boot_val()
    if (!is.null(b)) {
      key <- paste(b$Term, b$Type); m <- match(paste(r$Term, r$Type), key)
      r$lo <- b$CI_lower[m]; r$hi <- b$CI_upper[m]; ci_lab <- "bootstrap 95% CI"
    } else {
      r$lo <- r$Estimate - 1.96 * r$Std_Error; r$hi <- r$Estimate + 1.96 * r$Std_Error; ci_lab <- "± 1.96 SE"
    }
    r$Type <- factor(r$Type, levels = c("Linear", "Quadratic", "Correlational"),
                     labels = c("Directional (β)", "Quadratic (γ)", "Correlational (γij)"))
    r$Term <- factor(r$Term, levels = rev(unique(r$Term)))
    ggplot(r, aes(x = Estimate, y = Term)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, orientation = "y") +
      geom_point(size = 3, colour = "#2c7fb8") +
      facet_grid(Type ~ ., scales = "free_y", space = "free_y") +
      labs(x = paste0("Selection gradient (", ci_lab, ")"), y = NULL) + theme_bw() +
      theme(strip.text.y = element_text(angle = 0))
  })
  output$grad_plot <- renderPlot(grad_plot_obj())
  output$boot_note <- renderText({
    b <- boot_val()
    if (is.null(b)) "Error bars: ± 1.96 SE."
    else sprintf("Bootstrap: %d resamples, seed %d.", attr(b, "n_boot"), setup()$seed)
  })
  dl_name <- function(...) paste(c(file_stem(setup()$source), ...), collapse = "_")

  output$dl_report <- downloadHandler(
    filename = function() paste0(dl_name("selection_gradients"), ".csv"),
    content = function(f) {
      s <- setup(); r <- s$report; b <- boot_val()
      if (!is.null(b)) { m <- match(paste(r$Term, r$Type), paste(b$Term, b$Type)); r$Boot_SE <- b$Boot_SE[m]; r$CI_lower <- b$CI_lower[m]; r$CI_upper <- b$CI_upper[m] }
      utils::write.csv(as.data.frame(r), f, row.names = FALSE)
    })
  output$dl_gradplot <- downloadHandler(
    filename = function() paste0(dl_name("selection_gradients"), ".png"),
    content = function(f) ggsave(f, grad_plot_obj(), width = 7, height = 4.5, dpi = 200))

  # ---- Fitness functions tab ----
  # The trait in use is held in a reactiveVal: Run fits the default trait at
  # once, and the select's round trip with the same value does not refit it.
  uni_choice <- reactiveVal(NULL)
  observeEvent(input$uni_trait, {
    s <- setup()
    if (isTRUE(input$uni_trait %in% s$traits)) uni_choice(input$uni_trait)
  })
  uni_fit <- reactive({
    s <- setup(); tr <- uni_choice(); req(tr)
    with_seed(s$seed, suppressWarnings(suppressMessages(
      univariate_spline(s$prep, s$fit, tr, fitness_type = s$ftype, group = s$group_model,
                        k = input$spline_k, bs = s$spline_bs, smoothing = s$spline_sm,
                        bootstrap = TRUE, n_boot = 200, count_family = s$count_family))))
  })
  uni_plot_obj <- reactive(plot_univariate_fitness(uni_fit(), uni_choice(), classic_plot = input$classic) +
                             labs(x = sd_lab(uni_choice())))
  output$uni_plot <- renderPlot(uni_plot_obj())
  output$dl_uniplot <- downloadHandler(
    filename = function() paste0(dl_name("fitness_function", uni_choice()), ".png"),
    content = function(f) ggsave(f, uni_plot_obj(), width = 7, height = 5, dpi = 200))
  uni_land <- reactive({
    s <- setup(); u <- uni_fit()
    with_seed(s$seed, suppressWarnings(suppressMessages(capture.output(
      out <- adaptive_landscape(s$prep, u$model, uni_choice(), group_col = s$group_model,
                                grid_n = s$grid_n, simulation_n = s$sim_n, clamp = s$clamp)))))
    out
  })
  uni_land_plot_obj <- reactive({
    req(isTRUE(input$uni_land))
    # traits are standardised (within group), so the current mean is 0
    plot_adaptive_landscape(uni_land(), uni_choice(), show_actual_means = FALSE) +
      geom_vline(xintercept = 0, colour = "#e74c3c", alpha = 0.6) +
      annotate("text", x = 0, y = Inf, label = "current mean", vjust = 1.5, hjust = -0.1,
               size = 3.2, colour = "#e74c3c") +
      labs(x = sd_lab(uni_choice(), mean = TRUE))
  })
  output$uni_land_plot <- renderPlot(uni_land_plot_obj())
  output$dl_unilandplot <- downloadHandler(
    filename = function() paste0(dl_name("adaptive_landscape", uni_choice()), ".png"),
    content = function(f) ggsave(f, uni_land_plot_obj(), width = 7, height = 5, dpi = 200))

  # ---- Fitness surface and landscape ----
  # Same idea as uni_choice: the pair in use lives in a reactiveVal
  surf_choice <- reactiveVal(NULL)
  pick_surface <- function(s, x, y) {
    if (is.null(x) || !x %in% s$traits) x <- s$traits[1]
    if (is.null(y) || !y %in% s$traits || y == x) y <- setdiff(s$traits, x)[1]
    c(x, y)
  }
  observeEvent(c(input$surf_x, input$surf_y), {
    s <- setup(); surf_choice(pick_surface(s, input$surf_x, input$surf_y))
  })
  surf_traits <- reactive({
    tr <- surf_choice(); req(tr)
    validate(need(!is.na(tr[2]), "Pick two different traits for a surface."))
    tr
  })
  surfaces <- reactive({
    s <- setup(); tr <- surf_traits()
    if (s$surf_method == "tps" && !has_fields) validate("Thin-plate spline needs the 'fields' package; switch to GAM in Advanced settings.")
    surf <- suppressWarnings(suppressMessages(
      correlated_fitness_surface(s$prep, s$fit, tr, method = s$surf_method, grid_n = s$surf_grid, mask = !s$surf_full,
                                 too_far = s$surf_far, group = s$group, group_effect = s$within_group,
                                 k = s$surf_k, bs = s$surf_bs, smoothing = s$surf_sm, clamp = s$clamp,
                                 count_family = s$count_family)))
    land <- with_seed(s$seed, suppressWarnings(suppressMessages(capture.output(
      out <- adaptive_landscape(s$prep, surf$model, tr, group_col = s$group_model,
                                grid_n = s$grid_n, simulation_n = s$sim_n, clamp = s$clamp)))))
    list(surf = surf, land = out, traits = tr)
  })
  surf_plot_obj <- reactive({
    sf <- surfaces(); s <- setup()
    groups <- isTRUE(input$show_groups)
    lines <- !isFALSE(input$group_lines)
    alpha <- if (is.numeric(input$point_alpha)) input$point_alpha else 0.5
    fill_name <- if (s$ftype == "binary") "Survival" else "Fitness"
    # the thin-plate spline has no standard errors to draw
    unc <- if (s$surf_method == "tps" || is.null(input$surf_unc)) "none" else input$surf_unc
    p <- if (isTRUE(input$show_points)) {
      suppressMessages(plot_correlated_fitness_enhanced(sf$surf, sf$traits, original_data = s$prep, fitness_col = s$fit, bins = 12,
                                                        point_alpha = alpha, show_groups = groups, group_lines = lines,
                                                        uncertainty = unc, fill = fill_name))
    } else plot_correlated_fitness(sf$surf, sf$traits, bins = 12, show_groups = groups, group_lines = lines,
                                   uncertainty = unc, fill = fill_name)
    p <- suppressMessages(apply_theme(p, input$theme, binary = if (isTRUE(input$show_points)) s$ftype == "binary" else NULL,
                                      fill_name = fill_name, point_name = s$fit))
    p + labs(x = sd_lab(sf$traits[1]), y = sd_lab(sf$traits[2]))
  })
  output$surf_plot <- renderPlot(surf_plot_obj())
  output$dl_surfplot <- downloadHandler(
    filename = function() paste0(dl_name("fitness_surface", surf_traits()), ".png"),
    content = function(f) ggsave(f, surf_plot_obj(), width = 7.5, height = 6, dpi = 200))

  output$opt_txt <- renderText({
    sf <- surfaces(); opt <- sf$land$optimum; tr <- sf$traits; g <- sf$land$grid
    dx <- opt[[tr[1]]]; dy <- opt[[tr[2]]]
    on_edge <- vapply(tr, function(t) {
      v <- g[[t]]; opt[[t]] <= min(v) + 1e-9 || opt[[t]] >= max(v) - 1e-9
    }, logical(1))
    where <- sprintf("%s = %+.2f SD and %s = %+.2f SD from the current population mean", tr[1], dx, tr[2], dy)
    sup <- sf$land$support
    beyond <- if (is.null(sup)) "" else sprintf(" %.0f%% of the population simulated there lies outside the data.", 100 * sup$at_optimum)
    if (any(on_edge)) {
      sprintf("No interior optimum: mean fitness keeps rising towards %s (highest %.3f).%s", where, opt$.mean_fit, beyond)
    } else {
      sprintf("Optimum at %s, %.2f SD from the current mean (mean fitness %.3f).%s", where, sqrt(dx^2 + dy^2), opt$.mean_fit, beyond)
    }
  })
  land_plot_obj <- reactive({
    sf <- surfaces(); tr <- sf$traits
    p <- plot_adaptive_landscape(sf$land, tr, bins = 12,
                                 show_optimum = isTRUE(input$show_opt), show_actual_means = FALSE,
                                 show_support = isTRUE(input$show_support))
    p <- suppressMessages(apply_theme(p, input$land_theme, fill_name = "Mean fitness")) +
      labs(x = sd_lab(tr[1], mean = TRUE), y = sd_lab(tr[2], mean = TRUE))
    if (isTRUE(input$show_mean)) {
      # current population mean is (0, 0) in SD units
      p <- p + annotate("point", x = 0, y = 0, shape = 21, size = 3.5, fill = "#e74c3c", colour = "black") +
        annotate("text", x = 0, y = 0, label = "Current mean", vjust = -1.1, size = 3)
    }
    p
  })
  output$land_plot <- renderPlot(land_plot_obj())
  output$land3d_ui <- renderUI({
    if (has_plotly) plotly::plotlyOutput("land3d_plotly", height = "560px")
    else plotOutput("land3d", height = "560px")
  })
  output$land3d <- renderPlot({
    sf <- surfaces()
    if (!has_fields) { plot.new(); text(0.5, 0.5, "3D view needs the 'fields' or 'plotly' package"); return() }
    plot_adaptive_landscape_3d(sf$land, sf$traits, color_palette = theme_of(input$land_theme)$ramp(100))
  })
  # plotly is optional: without it the 3D view above falls back to fields
  if (has_plotly) output$land3d_plotly <- plotly::renderPlotly({
    sf <- surfaces(); g <- sf$land$grid; tr <- sf$traits; opt <- sf$land$optimum
    xu <- sort(unique(g[[tr[1]]])); yu <- sort(unique(g[[tr[2]]]))
    z <- matrix(NA_real_, nrow = length(yu), ncol = length(xu)) # plotly: rows follow y
    z[cbind(match(g[[tr[2]]], yu), match(g[[tr[1]]], xu))] <- g$.mean_fit
    pal <- theme_of(input$land_theme)$ramp(9)
    cs <- lapply(seq_along(pal), function(i) list((i - 1) / (length(pal) - 1), pal[i]))
    # current population mean is (0, 0) in SD units
    i0 <- which.min((g[[tr[1]]])^2 + (g[[tr[2]]])^2)
    p <- plotly::plot_ly()
    p <- plotly::add_surface(p, x = xu, y = yu, z = z, colorscale = cs, opacity = 0.95,
                             colorbar = list(title = "Mean fitness"),
                             hovertemplate = paste0(tr[1], ": %{x:.2f}<br>", tr[2], ": %{y:.2f}<br>mean fitness: %{z:.3f}<extra></extra>"))
    if (isTRUE(input$show_opt)) {
      p <- plotly::add_markers(p, x = opt[[tr[1]]], y = opt[[tr[2]]], z = opt$.mean_fit,
                               marker = list(size = 6, color = "#f1c40f", symbol = "diamond", line = list(color = "black", width = 1)),
                               name = "Optimum", hoverinfo = "name")
    }
    if (isTRUE(input$show_mean)) {
      p <- plotly::add_markers(p, x = g[[tr[1]]][i0], y = g[[tr[2]]][i0], z = g$.mean_fit[i0],
                               marker = list(size = 5, color = "#e74c3c", symbol = "circle", line = list(color = "black", width = 1)),
                               name = "Current population mean", hoverinfo = "name")
    }
    plotly::layout(p,
      title = list(text = "Adaptive landscape (drag to rotate)", font = list(size = 14)),
      showlegend = FALSE, margin = list(l = 0, r = 0, b = 0, t = 40),
      scene = list(
        xaxis = list(title = sd_lab(tr[1], mean = TRUE)), yaxis = list(title = sd_lab(tr[2], mean = TRUE)),
        zaxis = list(title = "Mean fitness"),
        camera = list(eye = list(x = 1.6, y = -1.6, z = 0.9))))
  })
  output$dl_landplot <- downloadHandler(
    filename = function() paste0(dl_name("adaptive_landscape", surfaces()$traits), ".png"),
    content = function(f) ggsave(f, land_plot_obj(), width = 7.5, height = 6, dpi = 200))

  # ---- Groups tab ----
  output$groups_ui <- renderUI({
    s <- setup()
    if (is.null(s$group)) return(div(class = "help-note", "Pick a group and tick 'Estimate per group as well'."))
    if (is.null(s$grouped)) return(div(class = "help-note", "Tick 'Estimate per group as well' and run again."))
    tagList(
      res_head(sprintf("Selection gradients by %s", s$group),
               dl_btn("dl_groups", "CSV", "Download the gradients by group as CSV")),
      tableOutput("group_table"),
      res_head("Gradient plot by group", dl_btn("dl_groupplot", "PNG", "Download the plot by group as PNG")),
      plotOutput("group_plot", height = "380px"),
      div(class = "help-note", "Each group standardised on its own.")
    )
  })
  output$group_table <- renderTable({
    s <- setup(); g <- s$grouped
    req(!is.null(g))
    data.frame(Group = as.character(g$Group), Term = g$Term, Type = g$Type,
               Estimate = fmt(g$Beta_Coefficient), SE = fmt(g$Standard_Error),
               p = paste(sub("^= ", "", fmt_p(g$P_Value)), stars(g$P_Value)),
               check.names = FALSE, stringsAsFactors = FALSE)
  }, align = "lllrrr")
  output$dl_groups <- downloadHandler(
    filename = function() paste0(dl_name("selection_gradients_by", setup()$group), ".csv"),
    content = function(f) {
      g <- setup()$grouped; req(!is.null(g))
      utils::write.csv(as.data.frame(g), f, row.names = FALSE)
    })
  group_plot_obj <- reactive({
    s <- setup(); g <- s$grouped
    req(!is.null(g))
    g <- g[g$Type != "Correlational" | length(s$traits) <= 3, ]
    g$Group <- factor(g$Group)
    g$Type <- factor(g$Type, levels = c("Linear", "Quadratic", "Correlational"),
                     labels = c("Directional (β)", "Quadratic (γ)", "Correlational (γij)"))
    ggplot(g, aes(x = Beta_Coefficient, y = Term, colour = Group)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_errorbar(aes(xmin = Beta_Coefficient - 1.96 * Standard_Error, xmax = Beta_Coefficient + 1.96 * Standard_Error),
                    width = 0.25, orientation = "y", position = position_dodge(width = 0.5)) +
      geom_point(size = 3, position = position_dodge(width = 0.5)) +
      facet_grid(Type ~ ., scales = "free_y", space = "free_y") +
      labs(x = "Selection gradient (± 1.96 SE)", y = NULL) + theme_bw() +
      theme(strip.text.y = element_text(angle = 0))
  })
  output$group_plot <- renderPlot(group_plot_obj())
  output$dl_groupplot <- downloadHandler(
    filename = function() paste0(dl_name("selection_gradients_by", setup()$group), ".png"),
    content = function(f) ggsave(f, group_plot_obj(), width = 7, height = 4.5, dpi = 200))

  # ---- Settings used (Data tab) ----
  code_lines <- reactive({
    # the settings of the results on screen: the data and traits of the last
    # run, the pair the surface was drawn for, the resamples drawn
    s <- setup(); b <- boot_val(); tr <- surf_choice()
    r_code(s, s$dataset, s$file_name,
           uni_trait = uni_choice() %||% s$traits[1],
           spline_k = input$spline_k %||% 10,
           surf_traits = if (length(s$traits) >= 2 && length(tr) == 2 && !anyNA(tr) && tr[1] != tr[2]) tr else character(),
           n_boot = if (!is.null(b)) attr(b, "n_drawn") else if (is.numeric(input$n_boot) && !is.na(input$n_boot)) round(input$n_boot) else 500,
           uncertainty = if (s$surf_method == "tps" || is.null(input$surf_unc)) "none" else input$surf_unc,
           canonical = isTRUE(input$canonical), uni_land = isTRUE(input$uni_land))
  })
  output$r_code <- renderText(paste(code_lines(), collapse = "\n"))
  output$dl_code <- downloadHandler(
    filename = function() paste0(dl_name("selection_analysis"), ".R"),
    content = function(f) writeLines(code_lines(), f))

  output$settings <- renderText({
    s <- setup(); b <- boot_val()
    # mgcv turns GCV.Cp into UBRE for survival and Poisson counts and into REML
    # for a negative binomial; the fitness function takes the app's fitness
    # type, the surface finds the type in the data
    crit <- function(sm, type) if (!identical(sm, "GCV.Cp")) sm else if (type == "binary") "UBRE" else
      if (type == "count") switch(s$count_family, poisson = "UBRE", nb = "REML", "GCV") else "GCV"
    surf_type <- suppressWarnings(detect_family(s$d[[s$fit]])$type)
    surf_crit <- crit(s$surf_sm, surf_type)
    paste(
      sprintf("File: %s", s$source),
      sprintf("Fitness: %s (%s, %s)", s$fit, s$ftype, s$ftype_how),
      sprintf("Traits: %s", paste(s$traits, collapse = ", ")),
      if (is.null(s$group)) "Group: none" else
        sprintf("Group: %s (%s%s)", s$group,
                if (s$within_group) "standardised and relativised within group, group as a fixed effect" else
                  "one population, standardised together; group marks means and peaks on the surface",
                if (s$per_group) "; selection also estimated per group" else ""),
      sprintf("Individuals with complete data: %d", nrow(s$d)),
      sprintf("Fitness function: %s basis, %s smoothing, k = %d%s", s$spline_bs,
              crit(s$spline_sm, s$ftype), input$spline_k,
              if (identical(s$ftype, "count")) paste0(", ", s$count_family, " family") else ""),
      sprintf("Fitness surface: %s, %d x %d grid, %s",
              if (s$surf_method == "tps") paste0("thin-plate spline", if (s$clamp) ", held within the fitness range" else "") else
                sprintf("GAM with %s basis, %s smoothing, k %s%s", s$surf_bs, surf_crit,
                        if (is.null(s$surf_k)) "from the data" else paste("=", s$surf_k),
                        if (surf_type == "count") paste0(", ", s$count_family, " family") else ""),
              s$surf_grid, s$surf_grid, paste0(if (s$surf_full) "drawn over the full grid" else "blank outside the data",
                if (!is.null(s$surf_far)) sprintf(", cells farther than %s of the range from any individual blank", s$surf_far) else "",
                if (!is.null(s$group)) paste0(", group means and peaks marked", if (s$within_group) ", group as a fixed effect" else ", one surface for all groups") else "")),
      sprintf("Adaptive landscape: %d x %d grid, %d simulated individuals per point", s$grid_n, s$grid_n, s$sim_n),
      if (is.null(b)) "Bootstrap: not run" else sprintf("Bootstrap: %d resamples", attr(b, "n_boot")),
      sprintf("Random seed: %d", s$seed),
      sprintf("lande %s, R %s.%s", as.character(utils::packageVersion("lande")), R.version$major, R.version$minor),
      sep = "\n")
  })
}

shinyApp(ui, server)
