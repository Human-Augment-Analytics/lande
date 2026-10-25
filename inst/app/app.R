# ---------------------------------------------------------------------------
# lande: interactive selection analysis
#
# Run locally:   shiny::runApp("app")
# Static export (GitHub Pages): shinylive::export("app", "docs"); serve docs/
#
# The core analysis needs only mgcv and MASS. The thin-plate spline surface,
# the 3D landscape and surface (fields) and Type III ANOVA (car) degrade to a
# message when those packages are missing.
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
has_zip <- requireNamespace("zip", quietly = TRUE)
# minutes without a click or keypress before a hosted session is dropped; 0
# turns the disconnect off
IDLE_MINUTES <- getOption("lande.idle_minutes", 10)

# Maintenance notice. On the lab server the updater writes a small JSON file
# (restart time, reason) before it restarts the app, and starts the workers
# with LANDE_NOTICE_FILE pointing at it. Unset (the GitHub Pages copy, a local
# run) means none of this is used. A missing, partly written or malformed file
# counts as no notice, and so does one more than 15 minutes past its time.
NOTICE_FILE <- Sys.getenv("LANDE_NOTICE_FILE")
read_notice <- function(path) {
  off <- list(active = FALSE)
  txt <- tryCatch(suppressWarnings(readLines(path, n = 50L, warn = FALSE)), error = function(e) NULL)
  if (!length(txt)) return(off)
  x <- tryCatch(jsonlite::fromJSON(paste(txt, collapse = "\n")), error = function(e) NULL)
  at <- suppressWarnings(as.numeric(if (is.list(x)) x$restart_at))
  if (length(at) != 1L || !is.finite(at)) return(off)
  left <- at - as.numeric(Sys.time())
  if (left < -900) return(off)
  when <- if (is.character(x$at) && length(x$at) == 1L && nzchar(x$at)) substr(x$at, 1, 20)
          else format(as.POSIXct(at, origin = "1970-01-01"), "%H:%M")
  list(active = TRUE, restart_at = at, seconds = round(left), at = when, update = identical(x$reason, "update"))
}

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

# Each heavy step is stored under a hash of its data, every argument it is
# called with and, for the random steps, the seed. A Run, a tab or a slider
# that comes back to settings already fitted then gets the stored result
# instead of a refit; as each random step sets its own seed, that is the
# result a refit would give. A bundled dataset's results are shared by all
# the sessions of this R process, so a second visitor on a server gets them
# too; an uploaded file's stay in its own session.
SHARED_STORE <- cachem::cache_mem(max_size = 200 * 1024^2)
new_store <- function() cachem::cache_mem(max_size = 100 * 1024^2)
memo <- function(store, key, expr) {
  k <- rlang::hash(key)
  value <- store$get(k)
  if (cachem::is.key_missing(value)) {
    value <- expr
    store$set(k, value)
  }
  value
}

# The fitness-function band's 200 refits run on worker processes when the
# machine has cores and memory to spare: one socket cluster per R process,
# shared by its sessions, started at the first Run and stopped with the app.
# The band is the same either way. options(lande.workers = 0) keeps the
# refits in this process, and a number of 2 or more sets the workers. The
# browser version (webR), machines with fewer than three cores and Linux
# with too little memory always refit here, and so, unless lande.workers
# asks for workers, do Windows and the hosting services.
BAND <- new.env()
# the cores a container's CPU quota allows (Linux), NA without a quota
cpu_quota <- function() {
  v2 <- tryCatch(scan("/sys/fs/cgroup/cpu.max", what = "", quiet = TRUE), error = function(e) NULL, warning = function(w) NULL)
  if (length(v2) == 2 && v2[1] != "max") return(as.numeric(v2[1]) / as.numeric(v2[2]))
  read1 <- function(f) tryCatch(as.numeric(readLines(f, 1, warn = FALSE)), error = function(e) NA, warning = function(w) NA)
  q <- read1("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")
  if (isTRUE(q > 0)) q / read1("/sys/fs/cgroup/cpu/cpu.cfs_period_us") else NA
}
# the memory left for workers on Linux, in bytes: what a container's limit
# (cgroup v2 or v1) leaves above what it uses, or what the machine has
# available, whichever is less; NA elsewhere or when neither can be read
mem_room <- function() {
  if (!identical(Sys.info()[["sysname"]], "Linux")) return(NA_real_)
  read1 <- function(f) tryCatch(as.numeric(readLines(f, 1, warn = FALSE)), error = function(e) NA, warning = function(w) NA)
  room <- NA_real_
  for (cg in list(c("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory.current"),
                  c("/sys/fs/cgroup/memory/memory.limit_in_bytes", "/sys/fs/cgroup/memory/memory.usage_in_bytes"))) {
    limit <- read1(cg[1])
    # no limit reads as "max" (v2) or a number near 2^63 (v1)
    if (isTRUE(limit < 2^60)) {
      room <- limit - max(0, read1(cg[2]), na.rm = TRUE)
      break
    }
  }
  info <- tryCatch(readLines("/proc/meminfo", warn = FALSE), error = function(e) character())
  avail <- suppressWarnings(as.numeric(sub("^MemAvailable:\\s*([0-9]+) kB.*$", "\\1", grep("^MemAvailable:", info, value = TRUE)[1])) * 1024)
  if (isTRUE(avail >= 0)) room <- min(room, avail, na.rm = TRUE)
  room
}
band_workers <- function() {
  if (grepl("emscripten", R.version$platform, fixed = TRUE)) return(0L)
  n <- getOption("lande.workers")
  if (!is.null(n)) return(if (isTRUE(n >= 2)) as.integer(n) else 0L)
  # a socket cluster's master listens on every network interface, which on
  # Windows brings up the firewall's prompt at the first Run
  if (.Platform$OS.type == "windows") return(0L)
  # a hosting service runs several app processes on one instance, each of
  # which would start its own workers, mostly in too little memory
  if (Sys.getenv("R_CONFIG_ACTIVE") %in% c("shinyapps", "rsconnect")) return(0L)
  cores <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
  quota <- suppressWarnings(cpu_quota())
  if (isTRUE(quota < cores)) cores <- floor(quota)
  if (!isTRUE(cores >= 3)) return(0L)
  n <- min(4, cores - 1)
  # a worker takes about 230 MB once it has refitted the band; each gets
  # 300 MB, after 1 GB for this process to grow into
  room <- suppressWarnings(mem_room())
  if (!is.na(room)) n <- min(n, floor((room - 1024^3) / (300 * 1024^2)))
  if (n < 2) 0L else as.integer(n)
}
start_band_cluster <- function(n) {
  cl <- parallel::makePSOCKcluster(n, rscript_args = "--vanilla", setup_timeout = 10)
  ok <- FALSE
  on.exit(if (!ok) parallel::stopCluster(cl))
  # the workers look where this process does; the function goes without its
  # environment, which would carry this app's along
  set_libs <- function(paths) invisible(.libPaths(paths))
  environment(set_libs) <- baseenv()
  parallel::clusterCall(cl, set_libs, .libPaths())
  # they refit with mgcv alone (lande is not loaded on them), which has to
  # be the version this process has
  theirs <- parallel::clusterEvalQ(cl, getNamespaceVersion(loadNamespace("mgcv")))
  if (!all(vapply(theirs, identical, logical(1), getNamespaceVersion(loadNamespace("mgcv"))))) {
    stop("the workers load another version of mgcv")
  }
  ok <- TRUE
  cl
}
# the cluster, started on first use; NULL to refit here
band_cluster <- function() {
  if (is.null(BAND$cl) && is.null(BAND$off)) {
    n <- band_workers()
    BAND$cl <- if (n >= 2) tryCatch(start_band_cluster(n), error = function(e) {
      try(message("lande app: no workers for the band (", conditionMessage(e), "); refitting it in this process"), silent = TRUE)
      NULL
    })
    # one try per process: without workers the band is refitted here
    BAND$off <- is.null(BAND$cl)
  }
  BAND$cl
}
stop_band_cluster <- function() {
  if (!is.null(BAND$cl)) try(parallel::stopCluster(BAND$cl), silent = TRUE)
  BAND$cl <- NULL
}
onStop(stop_band_cluster)

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

# read an uploaded table as spreadsheets write it: UTF-8 (with or without a
# byte-order mark), UTF-16, or Windows' own encoding when the text isn't valid
# UTF-8 (Excel on Windows); commas, or semicolons with decimal commas (Excel
# in much of Europe). The read.csv arguments used go with the data, for the
# R script
read_upload <- function(path) {
  size <- file.size(path)
  if (is.na(size) || size == 0) stop("the file is empty")
  lead <- readBin(path, "raw", 3)
  enc <- if (length(lead) >= 2 && (identical(lead[1:2], as.raw(c(0xff, 0xfe))) || identical(lead[1:2], as.raw(c(0xfe, 0xff))))) {
    "UTF-16"
  } else if (identical(lead, as.raw(c(0xef, 0xbb, 0xbf)))) {
    "UTF-8-BOM"
  } else if (all(validUTF8(suppressWarnings(readLines(path, n = 500, warn = FALSE))))) {
    "UTF-8"
  } else {
    "windows-1252"
  }
  con <- file(path, encoding = enc)
  first <- tryCatch(suppressWarnings(readLines(con, n = 20, warn = FALSE)), finally = close(con))
  first <- first[nzchar(trimws(first))]
  if (!length(first)) stop("the file has no rows")
  # the separator is whichever of comma, semicolon or tab the header has most of
  n_of <- function(ch) lengths(regmatches(first[1], gregexpr(ch, first[1], fixed = TRUE)))
  counts <- c(comma = n_of(","), semicolon = n_of(";"), tab = n_of("\t"))
  args <- switch(names(which.max(counts)),
                 semicolon = list(sep = ";", dec = ","),
                 tab = list(sep = "\t", dec = "."),
                 list(sep = ",", dec = "."))
  args$fileEncoding <- enc
  # a blank cell is missing, as a spreadsheet writes it, not a group named ""
  d <- suppressWarnings(do.call(utils::read.csv, c(list(path, na.strings = c("NA", "")), args)))
  if (ncol(d) < 2) stop("it has only one column")
  if (!nrow(d)) stop("it has a header but no rows")
  list(data = d, read_args = args)
}

# zip a folder: the zip package where there is one (the browser has no zip
# program), R's own zip otherwise
zip_folder <- function(zipfile, dir) {
  if (has_zip) return(invisible(zip::zip(zipfile, basename(dir), root = dirname(dir))))
  owd <- setwd(dirname(dir)); on.exit(setwd(owd))
  stopifnot(utils::zip(zipfile, basename(dir), flags = "-r9Xq") == 0)
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

# the study each bundled dataset comes from, shown under the Dataset picker and
# in the settings so it can be cited; Bumpus (1899) has no DOI
martin_2016 <- list(
  short = "Martin, C. H. (2016) Evolution 70, 1265-1282.",
  cite = paste("Martin, C. H. (2016). Context dependence in complex adaptive landscapes: frequency and",
               "trait-dependent selection surfaces within an adaptive radiation of Caribbean pupfishes.",
               "Evolution 70, 1265-1282."),
  links = c("Paper" = "https://doi.org/10.1111/evo.12932",
            "Data (Dryad)" = "https://doi.org/10.5061/dryad.n3mj3"))
DATA_SOURCES <- list(
  "Bumpus sparrows" = list(
    short = "Bumpus, H. C. (1899) Biological Lectures, Woods Hole, 209-226.",
    cite = paste("Bumpus, H. C. (1899). The elimination of the unfit as illustrated by the introduced sparrow,",
                 "Passer domesticus. Biological Lectures, Marine Biological Laboratory, Woods Hole, 209-226."),
    links = character()),
  "Crescent Pond pupfish" = martin_2016,
  "Little Lake pupfish" = martin_2016,
  "Finch community (five groups)" = list(
    short = "Beausoleil, M.-O. et al. (2023) Evolution 77, 2533-2546.",
    cite = paste("Beausoleil, M.-O. et al. (2023). The fitness landscape of a community of Darwin's finches.",
                 "Evolution 77, 2533-2546."),
    links = c("Paper" = "https://doi.org/10.1093/evolut/qpad160",
              "Data (Borealis, CC BY 4.0)" = "https://doi.org/10.5683/SP3/0YIWSE")))
numeric_cols <- function(df) names(df)[vapply(df, is.numeric, logical(1))]

# Limits for a shared server, where one R process serves several people and a
# long computation holds up all of them. A fit takes longer with more
# individuals and more terms (k traits give k linear and k(k+1)/2 quadratic and
# correlational terms, and a group of g levels in the model gives g - 1 more),
# so the traits (all nine of Bumpus's fit), individuals and groups are capped,
# and so is the bootstrap: individuals x terms x resamples is at most what the
# finch data (3428 individuals, two traits) take at 5000 resamples. The groups
# matter most for the fitness function, whose band refits the curve 200 times
# at each Run, and its time grows faster than the number of groups. The limits
# are on when the app runs on its own, as on a server; run_app() turns them off
# for a session on your own machine, and options(lande.limits = ) sets them
# either way.
LIMITS <- isTRUE(getOption("lande.limits", TRUE))
MAX_TRAITS <- 10
MAX_ROWS <- 5000
MAX_GROUPS <- 20
max_boot <- function(n, k, groups = 1) {
  if (!LIMITS) return(5000)
  terms <- k + k * (k + 1) / 2 + groups - 1
  max(2, min(5000, floor(3428 * 5 * 5000 / (n * terms))))
}

# --- formatting helpers -----------------------------------------------------
stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
}
fmt <- function(x, d = 3) ifelse(is.na(x), "", formatC(x, digits = d, format = "f"))
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "< 0.001", paste("=", formatC(p, digits = 3, format = "f"))))
pm <- function(est, se) paste0(fmt(est), " ± ", fmt(se))
sq <- function(t) paste0(t, "²")

# the classic jet rainbow of many published fitness surfaces, dark blue through
# cyan, yellow and red to dark red; harder to read for colour-blind viewers
jet <- grDevices::colorRampPalette(c("#00007F", "#0000FF", "#007FFF", "#00FFFF", "#7FFF7F",
                                     "#FFFF00", "#FF7F00", "#FF0000", "#7F0000"))
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
                                    points = c("#D55E00", "#0072B2")),
  "Rainbow (jet), white and black" = list(fill = function(...) ggplot2::discrete_scale("fill", palette = jet, ...),
                                          ramp = jet,
                                          points = c("white", "black"))
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
  # encodeString escapes quotes and backslashes, so any name stays valid R
  if (is.null(x)) "NULL" else if (is.character(x)) {
    if (length(x) == 1) encodeString(x, quote = '"') else sprintf("c(%s)", paste(encodeString(x, quote = '"'), collapse = ", "))
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
  # an upload is read with the separator and encoding the app found in it
  load <- if (dataset %in% names(DATA_CODE)) DATA_CODE[[dataset]] else {
    a <- s$read_args %||% list()
    extra <- c(if (identical(a$sep, ";")) 'sep = ";", dec = ","', if (identical(a$sep, "\t")) 'sep = "\\t"',
               if (!is.null(a$fileEncoding) && a$fileEncoding != "UTF-8") paste("fileEncoding =", r_value(a$fileEncoding)))
    sprintf('dat <- read.csv(%s, na.strings = c("NA", "")%s)', r_value(file_name %||% "your_file.csv"),
            paste0(", ", extra, collapse = ""))
  }
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
                           simulation_n = s$sim_n, clamp = r_value(s$clamp), method = r_value(s$land_method),
                           assign = "land1"),
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
             simulation_n = s$sim_n, clamp = r_value(s$clamp), method = r_value(s$land_method), assign = "land"),
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
# the site's pages, so the app links back to all of them
SITE <- "https://lande.thestroudlab.com/"
CASE_STUDIES <- c("Bumpus sparrows" = "articles/bumpus-sparrows.html",
                  "San Salvador pupfish" = "articles/pupfish.html",
                  "Darwin's finch community" = "articles/finch-community.html",
                  "Darwin's finch annual selection" = "articles/finches-year-by-year.html")
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
     .well{background:#F1F4F0;border-color:#DEE5DF;box-shadow:none}
     .nav-tabs{border-bottom-color:#DEE5DF}
     .nav-tabs>li>a{color:#356653}
     .nav-tabs>li>a:hover,.nav-tabs>li>a:focus{background:#EEF3EF;border-color:#EEF3EF #EEF3EF #DEE5DF;color:#2B5444}
     .nav-tabs>li.active>a,.nav-tabs>li.active>a:hover,.nav-tabs>li.active>a:focus{color:#24332D;background:#FFFFFF;border-color:#DEE5DF #DEE5DF transparent;box-shadow:inset 0 3px 0 #356653}
     .tab-content{background:#FFFFFF;border:1px solid #DEE5DF;border-top:none;border-radius:0 0 4px 4px;padding:0 12px 12px}
     .btn-primary,.btn-primary:focus{background:#356653;border-color:#2B5444}
     .btn-primary:hover,.btn-primary:active,.btn-primary.active{background:#2B5444;border-color:#22433A}
     .form-control,.selectize-input{border-color:#7D8F85}
     .form-control:focus,.selectize-input.focus{border-color:#356653;box-shadow:0 0 0 2px rgba(53,102,83,.25)}
     .selectize-dropdown .active{background:#E3ECE6;color:#24332D}
     .selectize-control.multi .selectize-input>div{background:#E3ECE6;color:#24332D}
     .irs--shiny .irs-bar{background:#5E8C78;border-color:#5E8C78}
     .irs--shiny .irs-single,.irs--shiny .irs-from,.irs--shiny .irs-to{background:#356653}
     input[type=checkbox],input[type=radio]{accent-color:#356653}
     .progress-bar{background-color:#356653}
     .help-note{color:#4F5D57;font-size:13px;margin-top:6px}
     .help-note.warn{color:#8a4b00}
     .interp{background:#EEF3EF;border-left:3px solid #356653;padding:10px 14px;margin:10px 0;font-size:14px}
     .interp.caveat{background:#FBF6EA;border-left-color:#B7862B}
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
     .app-header .app-links a{text-decoration:underline;text-decoration-color:rgba(53,102,83,.45);text-underline-offset:3px}
     .app-header .app-links a:hover,.app-header .app-links a:focus{text-decoration-color:currentColor}
     .app-links details.menu{position:relative}
     .app-links details.menu summary{margin:0;font-weight:400;color:#356653;list-style:none;text-decoration:underline;text-decoration-color:rgba(53,102,83,.45);text-underline-offset:3px}
     .app-links details.menu summary::-webkit-details-marker{display:none}
     .app-links details.menu summary::after{content:' \\25BE';text-decoration:none;display:inline-block;margin-left:2px}
     .app-links details.menu .menu-list{position:absolute;right:0;top:calc(100% + 6px);z-index:1000;background:#FFFFFF;border:1px solid #DEE5DF;border-radius:4px;box-shadow:0 4px 12px rgba(36,51,45,.12);padding:6px 0;min-width:240px}
     .app-links details.menu .menu-list a{display:block;padding:6px 14px;text-decoration:none;white-space:nowrap}
     .app-links details.menu .menu-list a:hover,.app-links details.menu .menu-list a:focus{background:#EEF3EF}
     .dl-all{margin-top:8px}
     .intro{background:#EEF3EF;border-left:3px solid #356653;padding:10px 14px;margin:0 0 12px;font-size:14px}
     .res-head{display:flex;align-items:flex-end;justify-content:space-between;flex-wrap:wrap;gap:8px;margin:18px 0 8px}
     .res-head h4.sec{margin:0}
     .res-head .res-dl{display:flex;align-items:center;flex-wrap:wrap;gap:6px}
     .res-head .res-dl .form-group{margin:0}
     .res-head .res-dl .radio-inline{padding-top:0}
     @media (max-width:767px){#land_plot{height:400px !important}}"
  ))),
  tags$head(tags$link(rel = "icon", type = "image/png", href = "favicon.png")),
  # the Case studies menu closes on a click outside it or on one of its links
  tags$head(tags$script(HTML(
    "document.addEventListener('click', function (e) {
       document.querySelectorAll('details.menu[open]').forEach(function (d) {
         if (!d.contains(e.target) || e.target.closest('a')) d.removeAttribute('open');
       });
     });"))),
  # drop the connection after IDLE_MINUTES without a click or keypress, so a tab
  # left open does not keep the server instance awake
  if (IDLE_MINUTES > 0) tags$head(tags$script(HTML(paste0("
    (function () {
      var idleMinutes = ", IDLE_MINUTES, ";
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
  ")))),
  # the maintenance notice: a bar across the top of the page while a restart is
  # planned. It reaches the page two ways: the server sends it when the notice
  # file changes, and the page itself fetches lande-notice.json (which nginx
  # serves from the same file) every 30 seconds, so it arrives even while this
  # worker is busy with a long run and cannot send anything. The page keeps the
  # minutes current itself. Only on a server that sets LANDE_NOTICE_FILE.
  if (nzchar(NOTICE_FILE)) tags$script(HTML("
    (function () {
      var bar = null, timer = null, cur = null;
      function hide() {
        cur = null;
        if (timer) { clearInterval(timer); timer = null; }
        if (bar) { bar.remove(); bar = null; }
      }
      function show() {
        var left = Math.round((cur.deadline - Date.now()) / 1000);
        if (left < -900) { hide(); return; }
        var mins = Math.floor(left / 60);
        var text = left > 0
          ? (cur.update ? 'lande will restart for an update' : 'lande will restart') + ' at ' + cur.at + ' (in ' +
            (mins < 1 ? 'less than a minute' : mins === 1 ? '1 minute' : mins + ' minutes') + '). ' +
            'Download any results you need before then; you will need to reload the page afterwards.'
          : (cur.update ? 'lande is restarting for an update now.' : 'lande is restarting now.') +
            ' Reload the page in a minute or two.';
        if (bar.textContent !== text) bar.textContent = text;
      }
      // n: restart_at (epoch), seconds left by the server's clock (below 0 once
      // the time has passed), at (the time to show), update; or null
      function set(n) {
        if (!n) { hide(); return; }
        if (cur && cur.restart_at === n.restart_at && cur.update === n.update) return;
        cur = { restart_at: n.restart_at, at: n.at, update: n.update, deadline: Date.now() + 1000 * n.seconds };
        if (!bar) {
          bar = document.createElement('div');
          bar.id = 'lande-notice';
          bar.setAttribute('role', 'status');
          bar.style.cssText = 'position:sticky;top:0;z-index:2000;background:#fff3cd;color:#5c4400;' +
            'border-bottom:1px solid #e3c766;padding:10px 20px;font-size:15px;font-weight:600';
          document.body.insertBefore(bar, document.body.firstChild);
        }
        show();
        if (!timer) timer = setInterval(function () { if (cur) show(); }, 5000);
      }
      Shiny.addCustomMessageHandler('lande_notice', function (n) { set(n && n.active ? n : null); });
      function look() {
        fetch('lande-notice.json', { cache: 'no-store' }).then(function (r) {
          // not 200: nothing serves the file here (no nginx in front), so leave it to the server
          if (r.status !== 200) return;
          var now = Date.parse(r.headers.get('Date') || '');
          if (isNaN(now)) now = Date.now();
          return r.json().then(function (x) {
            var at = x && typeof x.restart_at === 'number' ? x.restart_at : NaN;
            if (!isFinite(at) || at - now / 1000 < -900) { set(null); return; }
            set({ restart_at: at, seconds: Math.round(at - now / 1000), update: x.reason === 'update',
                  at: typeof x.at === 'string' && x.at ? x.at.slice(0, 20) : new Date(at * 1000).toTimeString().slice(0, 5) });
          });
        }).catch(function () {});
      }
      setTimeout(look, 3000);
      setInterval(look, 30000);
    })();
  ")),
  div(class = "app-header",
      tags$img(src = "logo.png", alt = ""),
      div(tags$h1("lande"), div(class = "app-sub", "Evolutionary selection analysis")),
      div(class = "app-links",
          doc_link("Home", SITE),
          doc_link("Guide", paste0(SITE, "articles/evolutionary-selection-analysis.html")),
          tags$details(class = "menu",
                       tags$summary("Case studies"),
                       div(class = "menu-list", lapply(names(CASE_STUDIES), function(n)
                         doc_link(n, paste0(SITE, CASE_STUDIES[[n]]))))),
          doc_link("Reference", paste0(SITE, "reference/")),
          doc_link("Report a problem", "https://github.com/Human-Augment-Analytics/lande/issues"))),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Data"),
      selectInput("dataset", "Dataset", DATASETS),
      uiOutput("data_source"),
      conditionalPanel("input.dataset == 'Upload CSV...'",
        fileInput("file", "CSV file", accept = c(".csv", ".tsv", ".txt"))),
      h4("Design"),
      selectInput("fitness", "Fitness", NULL),
      div(class = "help-note", textOutput("fitness_hint")),
      selectizeInput("traits", "Traits", NULL, multiple = TRUE, options = if (LIMITS) list(maxItems = MAX_TRAITS)),
      selectInput("group", "Group (optional)", "(none)"),
      checkboxInput("within_group", "Standardise within group, group as a fixed effect", TRUE),
      div(class = "help-note",
          "Untick to treat all groups as one population; the group then only marks means and peaks on the surface, as for several species on one surface."),
      checkboxInput("per_group", "Estimate per group as well", FALSE),
      actionButton("run", "Run analysis", class = "btn-primary", width = "100%"),
      conditionalPanel("input.run > 0", div(class = "dl-all",
        downloadButton("dl_all", "Download all (.zip)", class = "btn-block",
                       title = "Download the data, the R script, the tables and the figures as one zip",
                       `aria-label` = "Download the data, the R script, the tables and the figures as one zip"))),
      tags$details(
        tags$summary("Advanced settings"),
        selectInput("ftype", "Fitness type", c("detect automatically" = "auto", "binary (0/1)" = "binary", "count" = "count", "continuous" = "continuous")),
        radioButtons("surf_method", "Fitness-surface method",
                     c("GAM (mgcv)" = "gam", "Thin-plate spline (fields)" = "tps"), inline = TRUE),
        sliderInput("surf_grid", "Surface grid", 20, 80, 50, step = 10),
        sliderInput("grid_n", "Landscape grid", 15, 60, 40, step = 5),
        selectInput("land_method", "Landscape method",
                    c("Lattice (fast, same as simulate)" = "lattice", "Gaussian convolution" = "convolve",
                      "Simulate (original)" = "simulate"), selected = "lattice"),
        div(class = "help-note",
            "Gaussian convolution draws no individuals: the simulated individuals below do not apply to it, and its landscape does not change with the seed."),
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
          div(class = "help-note", "S total selection; β directional; γ quadratic curvature (negative is consistent with stabilising selection, positive with disruptive); γij correlational. * p < 0.05, ** < 0.01, *** < 0.001."),
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
            column(3, conditionalPanel("input.surf_view == '2d'", br(),
              checkboxInput("show_points", "Show individuals", TRUE))),
            column(3, br(), checkboxInput("show_groups", "Group means and peaks", TRUE))),
          fluidRow(
            column(5, conditionalPanel("input.surf_view == '2d' && input.show_points",
              sliderInput("point_alpha", "Opacity of the individuals", 0.05, 1, 0.5, step = 0.05))),
            column(4, selectInput("theme", "Colours", names(THEMES))),
            column(3, conditionalPanel("input.surf_view == '2d'",
              selectInput("surf_unc", "Uncertainty (GAM)",
                          c("none" = "none", "standard error lines" = "se", "lower, fit and upper" = "band"))))),
          res_head("Fitness surface",
                   radioButtons("surf_view", tags$span(class = "sr-only", "Surface view"),
                                c("2D" = "2d", "3D" = "3d"), inline = TRUE),
                   conditionalPanel("input.surf_view == '2d'",
                     dl_btn("dl_surfplot", "PNG", "Download the fitness surface as PNG"))),
          conditionalPanel("input.surf_view == '2d'",
            plotOutput("surf_plot", height = "500px"),
            div(class = "help-note", "A group's peak is the highest point of the surface within that group's own range: filled when it is a peak of the surface, open when the surface keeps rising past the group's range or the edge of the data.")),
          conditionalPanel("input.surf_view == '3d'",
            uiOutput("surf3d_ui"),
            div(class = "help-note", "The fitted fitness where there are data. Drag to rotate. Diamonds: interior peaks of the surface; with groups shown, circles mark the group means and squares their highest points, filled when that point is a peak of the surface and open when it is not."))),
        tabPanel("Adaptive landscape",
          br(), uiOutput("opt_box"),
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

  # Keep a numeric setting within the range the sidebar offers. The browser
  # holds a slider or numericInput to its bounds, but the value comes over the
  # WebSocket and a hand-made message need not respect them, so every setting
  # that sets the size of a run is rounded and clamped before it is used.
  clamp_num <- function(x, lo, hi, default) {
    x <- suppressWarnings(round(as.numeric(x)))[1]
    if (length(x) == 0L || is.na(x)) return(default)
    min(max(x, lo), hi)
  }

  # mgcv and MASS would load with the first fit; loading them once the page
  # is drawn takes that time out of the first Run
  session$onFlushed(function() for (p in c("mgcv", "MASS")) try(loadNamespace(p), silent = TRUE))

  # the study the chosen dataset comes from, with links to the paper and the data
  output$data_source <- renderUI({
    src <- DATA_SOURCES[[input$dataset]]
    if (is.null(src)) return(NULL)
    links <- lapply(seq_along(src$links), function(i)
      tagList(if (i > 1) " | ", doc_link(names(src$links)[i], src$links[[i]])))
    # the short form here; the full citation on hover and in the settings
    div(class = "help-note", title = src$cite, tags$b("Source: "), src$short,
        if (length(links)) tagList(" ", links))
  })

  # results stored for this session's uploads (bundled datasets use SHARED_STORE)
  session_store <- new_store()
  # every heavy step goes through stored(); fits_made lists the steps this
  # session fitted rather than took from a store, for the tests
  fits_made <- character()
  stored <- function(store, key, expr) memo(store, key, {
    value <- expr
    fits_made <<- c(fits_made, key[[1]])
    value
  })

  current <- reactive({
    if (input$dataset == "Upload CSV...") {
      req(input$file)
      # Read only a file Shiny itself received. Shiny puts each upload in its own
      # randomly named folder inside this R process's temp folder; a hand-made
      # WebSocket message could otherwise set datapath to any path the process
      # can read, so anything outside that temp folder (or missing) is refused.
      p <- tryCatch(normalizePath(input$file$datapath[1], mustWork = TRUE), error = function(e) NULL)
      req(p, startsWith(p, paste0(normalizePath(tempdir()), .Platform$file.sep)))
      got <- tryCatch(read_upload(p), error = function(e) e)
      validate(need(!inherits(got, "error"), sprintf(
        "Couldn't read %s as a CSV file: %s. Save it from your spreadsheet as CSV and upload it again.",
        input$file$name, if (inherits(got, "error")) conditionMessage(got) else "")))
      list(data = got$data, fitness = NULL, traits = NULL, group = NULL,
           source = input$file$name, read_args = got$read_args, file_path = p)
    } else c(load_dataset(input$dataset), source = paste(input$dataset, "(bundled)"))
  })

  # Defaults for an uploaded file: a survival-like or 0/1 column as fitness, a
  # year or site column as group, and measured traits (not ids or years) as traits
  observeEvent(current(), {
    # a file that couldn't be read says so in the fitness hint; nothing to set here
    d <- tryCatch(current(), error = function(e) NULL)
    req(d)
    nm <- names(d$data); nums <- numeric_cols(d$data)
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
      # no more groups than the app takes (30 without the limits), counting
      # blanks as one more, as the package does
      g <- g[vapply(g, function(v) n_unique(v) >= 2 && length(unique(d$data[[v]])) <= (if (LIMITS) MAX_GROUPS else 30), logical(1))]
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
    # one check at a time, so the message names the first problem
    validate(
      need(length(traits) >= 1, "Pick at least one trait."),
      need(!LIMITS || length(traits) <= MAX_TRAITS, sprintf("Pick at most %d traits.", MAX_TRAITS)),
      need(!fit %in% traits, "Fitness column can't also be a trait."),
      need(all(c(fit, traits, grp) %in% names(d)), "Columns not found.")
    )
    validate(
      need(is.numeric(d[[fit]]) || is.logical(d[[fit]]), "Fitness column must be numeric (0/1 survival or a numeric measure)."),
      need(all(vapply(d[traits], is.numeric, logical(1))), "Trait columns must be numeric.")
    )
    validate(need(sum(stats::complete.cases(d[, c(fit, traits)])) >= 20, "Need at least 20 individuals with complete data."))
    validate(need(all(vapply(d[traits], function(x) length(unique(x[!is.na(x)])) > 2, logical(1))),
                  "Each trait needs more than two distinct values; a 0/1 column cannot be a trait."))
    # rows need fitness and every trait; the package keeps unlabelled rows as one more group
    d <- d[stats::complete.cases(d[, c(fit, traits)]), ]
    validate(need(nrow(d) >= 20, "Need at least 20 individuals with complete data."),
             need(!LIMITS || nrow(d) <= MAX_ROWS, sprintf("This app takes at most %d individuals with complete data; for more, use the lande package in R, where run_app() has no limits.", MAX_ROWS)))
    # a group is a column in the models (with the switch on), and each group is
    # also fitted, marked and listed on its own, so a column of many values (an
    # id, a measurement) picked as the group would tie the process up for hours;
    # blanks are one more group
    n_groups <- if (is.null(grp)) 1L else nlevels(addNA(factor(d[[grp]]), ifany = TRUE))
    validate(need(!LIMITS || n_groups <= MAX_GROUPS,
                  sprintf("This app takes at most %d groups, and %s has %d among the individuals with complete data. Pick another group column, or (none).",
                          MAX_GROUPS, grp, n_groups)))
    resolved <- if (input$ftype == "auto") detect_family(d[[fit]])$type else input$ftype
    if (!resolved %in% c("binary", "count")) resolved <- "continuous"
    # with the switch off the groups are one population: pooled standardisation,
    # no group term, the group only marks means and peaks on the surface
    within <- !is.null(grp) && !isFALSE(input$within_group)
    grp_model <- if (within) grp else NULL
    prep <- suppressWarnings(suppressMessages(
      prepare_selection_data(d, fit, traits, standardize = TRUE, group = grp_model, add_relative = TRUE, na_action = "none")))
    # the steps' results are stored under the hash of the data they are given
    store <- if (input$dataset == "Upload CSV...") session_store else SHARED_STORE
    d_key <- rlang::hash(d)
    args <- list(fit, traits, fitness_type = resolved, standardize = TRUE, group = grp_model)
    fitted <- stored(store, list("report", d_key, args), with_warnings(do.call(selection_report, c(list(d), args))))
    report <- fitted$value
    warned <- fitted$warnings
    grouped <- NULL
    if (!is.null(grp) && isTRUE(input$per_group)) {
      args <- list(fit, traits, fitness_type = resolved, standardize = TRUE, group = grp, return_grouped = TRUE)
      # a failure is not stored, so another Run (or session) tries again
      per <- tryCatch(stored(store, list("per_group", d_key, args), with_warnings(
        do.call(selection_coefficients, c(list(d), args)))),
        error = function(e) list(value = NULL, warnings = paste("Per-group fits failed:", conditionMessage(e))))
      grouped <- per$value
      if (length(per$warnings)) warned <- c(warned, paste("Per group:", per$warnings))
    }
    keep_or <- function(x, default) if (isTRUE(x %in% traits)) x else default
    updateSelectInput(session, "uni_trait", choices = traits, selected = keep_or(input$uni_trait, traits[1]))
    updateSelectInput(session, "surf_x", choices = traits, selected = keep_or(input$surf_x, traits[1]))
    updateSelectInput(session, "surf_y", choices = traits, selected = keep_or(input$surf_y, traits[min(2, length(traits))]))
    list(d = d, prep = prep, fit = fit, traits = traits, group = grp, n_groups = n_groups, ftype = resolved,
         report = report, grouped = grouped, warnings = unique(warned), source = current()$source,
         # what was run, for the script, even if the inputs change afterwards
         dataset = input$dataset, file_name = if (!is.null(input$file)) input$file$name,
         read_args = current()$read_args,
         file_path = current()$file_path,
         ftype_how = if (input$ftype == "auto") "detected" else "set by hand",
         per_group = isTRUE(input$per_group),
         seed = if (is.numeric(input$seed) && !is.na(input$seed)) round(input$seed) else 1,
         surf_method = input$surf_method, surf_grid = clamp_num(input$surf_grid, 20, 80, 50), surf_full = isTRUE(input$surf_full),
         surf_k = if (is.numeric(input$surf_k) && !is.na(input$surf_k)) clamp_num(input$surf_k, 5, 60, 5L) else NULL,
         surf_far = if (is.numeric(input$surf_far) && !is.na(input$surf_far) && input$surf_far > 0) input$surf_far else NULL,
         clamp = !isFALSE(input$clamp),
         within_group = within, group_model = grp_model,
         surf_bs = input$surf_bs %||% "tp", surf_sm = input$surf_sm %||% "REML",
         spline_bs = input$spline_bs %||% "cr", spline_sm = input$spline_sm %||% "GCV.Cp",
         count_family = input$count_family %||% "poisson",
         grid_n = clamp_num(input$grid_n, 15, 60, 40), sim_n = clamp_num(input$sim_n, 100, 1000, 300),
         land_method = input$land_method %||% "lattice",
         store = store, d_key = d_key, prep_key = rlang::hash(prep))
  })

  boot_val <- reactiveVal(NULL)
  # Run fits the gradients, then sends the data and gradients to the page at
  # once; the fitness function, surface and landscape are fitted straight
  # after, with the progress bar still up, so the tabs only draw afterwards
  observeEvent(input$run, {
    progress <- Progress$new(session)
    progress$set(message = "Running the analysis", value = 0.05, detail = "gradients")
    s <- try(setup(), silent = TRUE)
    if (inherits(s, "try-error")) {
      progress$close()
      req(FALSE)
    }
    uni_choice(if (isTRUE(input$uni_trait %in% s$traits)) input$uni_trait else s$traits[1])
    surf_choice(pick_surface(s, input$surf_x, input$surf_y))
    # the k the slider is on, even if it moved less than 600 ms ago
    take_k(input$spline_k)
    session$onFlushed(function() {
      on.exit(try(progress$close(), silent = TRUE))
      withReactiveDomain(session, isolate({
        progress$set(value = 0.2, detail = "fitness function")
        try(uni_fit(), silent = TRUE)
        if (isTRUE(input$uni_land)) try(uni_land(), silent = TRUE)
        progress$set(value = 0.6, detail = "fitness surface and landscape")
        if (length(s$traits) >= 2) try(landscape(), silent = TRUE)
        progress$set(value = 1, detail = "done")
      }))
    })
  }, priority = 100)
  observeEvent(input$run, boot_val(NULL))
  observeEvent(input$run_boot, {
    s <- setup()
    validate(need(is.numeric(input$n_boot) && !is.na(input$n_boot) && input$n_boot >= 2,
                  "Resamples must be a whole number of at least 2."))
    asked <- clamp_num(input$n_boot, 2, 5000, 500)
    # the group's columns count where the bootstrap fits them, in the survival
    # and count models; continuous fitness is only standardised within groups
    groups <- if (!is.null(s$group_model) && s$ftype != "continuous") s$n_groups else 1
    n_boot <- min(asked, max_boot(nrow(s$d), length(s$traits), groups))
    args <- list(s$fit, s$traits, fitness_type = s$ftype, standardize = TRUE, group = s$group_model,
                 n_boot = n_boot)
    withProgress(message = "Bootstrapping...", {
      b <- stored(s$store, list("bootstrap", s$d_key, args, s$seed), with_seed(s$seed, suppressWarnings(suppressMessages(
        do.call(bootstrap_selection, c(list(s$d), args))))))
    })
    # the number drawn, for the script; attr(b, "n_boot") counts the usable ones
    attr(b, "n_drawn") <- n_boot
    attr(b, "n_asked") <- asked
    attr(b, "groups") <- groups
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
    args <- list(s$fit, s$traits, fitness_type = s$ftype, group = s$group_model)
    stored(s$store, list("assumptions", s$d_key, args),
           suppressWarnings(do.call(check_selection_assumptions, c(list(s$d), args))))
  })
  output$assump_table <- renderTable({
    chk <- assump_checks()
    data.frame(
      Check = chk$check,
      Statistic = ifelse(is.na(chk$statistic), "", formatC(chk$statistic, digits = 3, format = "g")),
      p = ifelse(is.na(chk$p_value), "", ifelse(chk$p_value < 0.001, "< 0.001", formatC(chk$p_value, digits = 3, format = "f"))),
      Note = chk$note, check.names = FALSE, stringsAsFactors = FALSE)
  }, align = "lrrl")
  write_assump <- function(f) {
    chk <- assump_checks()
    utils::write.csv(data.frame(check = chk$check, statistic = chk$statistic, p_value = chk$p_value,
                                note = chk$note), f, row.names = FALSE)
  }
  output$dl_assump <- downloadHandler(
    filename = function() paste0(dl_name("assumption_checks"), ".csv"),
    content = write_assump)
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
    args <- list(s$fit, s$traits, fitness_type = s$ftype, standardize = TRUE, group = s$group_model)
    # a failure is not stored, so another Run (or session) tries again
    tryCatch(stored(s$store, list("canonical", s$d_key, args, s$seed), with_seed(s$seed, suppressWarnings(suppressMessages(
      do.call(canonical_analysis, c(list(s$d), args)))))),
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
      theme(strip.text.y = element_text(angle = 0), strip.background = element_rect(fill = "#E6EBE7"))
  })
  output$grad_plot <- renderPlot(grad_plot_obj())
  output$boot_note <- renderText({
    b <- boot_val()
    if (is.null(b)) "Error bars: ± 1.96 SE."
    else paste0(sprintf("Bootstrap: %d resamples, seed %d.", attr(b, "n_boot"), setup()$seed),
                if (attr(b, "n_asked") > attr(b, "n_drawn"))
                  sprintf(" %d were asked for; with this many %s the app runs at most %d.", attr(b, "n_asked"),
                          if (attr(b, "groups") > 1) "individuals, traits and groups" else "individuals and traits", attr(b, "n_drawn")))
  })
  dl_name <- function(...) paste(c(file_stem(setup()$source), ...), collapse = "_")

  write_report <- function(f) {
    s <- setup(); r <- s$report; b <- boot_val()
    if (!is.null(b)) { m <- match(paste(r$Term, r$Type), paste(b$Term, b$Type)); r$Boot_SE <- b$Boot_SE[m]; r$CI_lower <- b$CI_lower[m]; r$CI_upper <- b$CI_upper[m] }
    utils::write.csv(as.data.frame(r), f, row.names = FALSE)
  }
  output$dl_report <- downloadHandler(
    filename = function() paste0(dl_name("selection_gradients"), ".csv"),
    content = write_report)
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
  # the slider sends each value it rests on while dragged; the band follows
  # it once it has stayed put for 600 ms, and a Run takes the value it is on
  # (the Run observer calls take_k too). Each take marks the band as changed,
  # even at the same k, as the slider itself does, so the fitness function is
  # redrawn for every move the band has not yet taken but not again for one a
  # Run already took; the plot's jitter then moves the session's random
  # numbers as it would without the wait
  k_taken <- reactiveVal(list(k = NULL, n = 0))
  band_k <- reactive(k_taken()$k)
  k_moves <- 0
  k_used <- 0
  take_k <- function(k) {
    k_used <<- k_moves
    k_taken(list(k = clamp_num(k, 3, 12, 5), n = isolate(k_taken()$n) + 1))
  }
  observeEvent(input$spline_k, k_moves <<- k_moves + 1, priority = 200)
  k_rested <- debounce(reactive(input$spline_k), 600)
  observeEvent(k_rested(), if (k_moves > k_used) take_k(k_rested()))
  band_key <- reactive({
    s <- setup(); tr <- uni_choice(); req(tr)
    list("band", s$prep_key,
         # count_family only reaches a count fit, so other fits share one entry
         args = list(s$fit, tr, fitness_type = s$ftype, group = s$group_model, k = band_k(), bs = s$spline_bs,
                     smoothing = s$spline_sm, bootstrap = TRUE, n_boot = 200,
                     count_family = if (identical(s$ftype, "count")) s$count_family else "poisson"),
         s$seed)
  })
  uni_fit <- reactive({
    s <- setup(); key <- band_key()
    stored(s$store, key, {
      # the refits go to the workers where there are any; the band is the same
      band <- function(cl) with_seed(s$seed, suppressWarnings(suppressMessages(
        do.call(univariate_spline, c(list(s$prep), key$args, list(cl = cl))))))
      cl <- band_cluster()
      if (is.null(cl)) band(NULL) else tryCatch(band(cl), error = function(e) {
        # refitted here; if that works, the workers were the problem, so they go
        out <- band(NULL)
        stop_band_cluster()
        BAND$off <- TRUE
        out
      })
    })
  })
  uni_plot_obj <- reactive(plot_univariate_fitness(uni_fit(), uni_choice(), classic_plot = input$classic) +
                             labs(x = sd_lab(uni_choice())))
  output$uni_plot <- renderPlot(uni_plot_obj())
  output$dl_uniplot <- downloadHandler(
    filename = function() paste0(dl_name("fitness_function", uni_choice()), ".png"),
    content = function(f) ggsave(f, uni_plot_obj(), width = 7, height = 5, dpi = 200))
  uni_land <- reactive({
    s <- setup(); u <- uni_fit()
    args <- list(uni_choice(), group_col = s$group_model, grid_n = s$grid_n, simulation_n = s$sim_n,
                 clamp = s$clamp, method = s$land_method)
    stored(s$store, list("uni_land", band_key(), args, s$seed), {
      with_seed(s$seed, suppressWarnings(suppressMessages(capture.output(
        out <- do.call(adaptive_landscape, c(list(s$prep, u$model), args))))))
      out
    })
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
  # the surface and the landscape are fitted apart, so the surface tab does
  # not wait for the landscape and a landscape setting does not refit the surface
  surface_key <- reactive({
    s <- setup(); tr <- surf_traits()
    # count_family only reaches a surface fitted to counts, so other surfaces
    # share one entry; the surface finds counts by this same rule (fitness on
    # the rows with both traits), and if it can't be applied the family stays
    counts <- tryCatch({
      y <- as.numeric(s$prep[[s$fit]])
      keep <- stats::complete.cases(y, as.numeric(s$prep[[tr[1]]]), as.numeric(s$prep[[tr[2]]]))
      identical(suppressWarnings(detect_family(y[keep]))$type, "count")
    }, error = function(e) TRUE)
    list("surface", s$prep_key,
         args = list(s$fit, tr, method = s$surf_method, grid_n = s$surf_grid, mask = !s$surf_full,
                     too_far = s$surf_far, group = s$group, group_effect = s$within_group,
                     k = s$surf_k, bs = s$surf_bs, smoothing = s$surf_sm, clamp = s$clamp,
                     count_family = if (counts) s$count_family else "poisson"))
  })
  surface_fit <- reactive({
    s <- setup(); key <- surface_key()
    if (s$surf_method == "tps" && !has_fields) validate("Thin-plate spline needs the 'fields' package; switch to GAM in Advanced settings.")
    stored(s$store, key, suppressWarnings(suppressMessages(
      do.call(correlated_fitness_surface, c(list(s$prep), key$args)))))
  })
  landscape <- reactive({
    s <- setup(); surf <- surface_fit()
    args <- list(surf_traits(), group_col = s$group_model, grid_n = s$grid_n, simulation_n = s$sim_n,
                 clamp = s$clamp, method = s$land_method)
    stored(s$store, list("landscape", surface_key(), args, s$seed), {
      with_seed(s$seed, suppressWarnings(suppressMessages(capture.output(
        out <- do.call(adaptive_landscape, c(list(s$prep, surf$model), args))))))
      out
    })
  })
  surfaces <- reactive(list(surf = surface_fit(), land = landscape(), traits = surf_traits()))
  surf_plot_obj <- reactive({
    sf <- list(surf = surface_fit(), traits = surf_traits()); s <- setup()
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
  output$surf3d_ui <- renderUI({
    if (has_plotly) plotly::plotlyOutput("surf3d_plotly", height = "500px")
    else plotOutput("surf3d", height = "500px")
  })
  output$surf3d <- renderPlot({
    surf <- surface_fit(); tr <- surf_traits()
    if (!has_fields) { plot.new(); text(0.5, 0.5, "3D view needs the 'fields' or 'plotly' package"); return() }
    plot_correlated_fitness_3d(surf, tr, color_palette = theme_of(input$theme)$ramp(100))
  })
  # plotly is optional: without it the 3D view above falls back to fields
  if (has_plotly) output$surf3d_plotly <- plotly::renderPlotly({
    surf <- surface_fit(); tr <- surf_traits(); s <- setup(); g <- surf$grid
    xu <- sort(unique(g[[tr[1]]])); yu <- sort(unique(g[[tr[2]]]))
    # masked cells stay NA and are left as gaps
    z <- matrix(NA_real_, nrow = length(yu), ncol = length(xu)) # plotly: rows follow y
    z[cbind(match(g[[tr[2]]], yu), match(g[[tr[1]]], xu))] <- g$.fit
    pal <- theme_of(input$theme)$ramp(9)
    cs <- lapply(seq_along(pal), function(i) list((i - 1) / (length(pal) - 1), pal[i]))
    z_name <- if (s$ftype == "binary") "Survival" else "Fitness"
    # thin dark grey lines at every grid value outline each cell, as in the drape
    # view, faint enough that the colours stay light
    cell_lines <- function(v) list(show = TRUE, start = v[1], end = v[length(v)] + diff(v[1:2]) / 2,
                                   size = diff(v[1:2]), color = "rgba(77,77,77,0.5)", width = 1)
    p <- plotly::plot_ly()
    p <- plotly::add_surface(p, x = xu, y = yu, z = z, colorscale = cs, opacity = 0.95,
                             contours = list(x = cell_lines(xu), y = cell_lines(yu)),
                             # soft light, as the drape view's light shading: the colours
                             # stay near the colour bar and the slopes only dim a little
                             lighting = list(ambient = 0.9, diffuse = 0.35, specular = 0.02,
                                             roughness = 0.9, fresnel = 0.05),
                             colorbar = list(title = z_name),
                             hovertemplate = paste0(tr[1], ": %{x:.2f}<br>", tr[2], ": %{y:.2f}<br>", tolower(z_name), ": %{z:.3f}<extra></extra>"))
    pk <- surf$peaks
    pk <- if (is.null(pk)) NULL else pk[pk$interior %in% TRUE, , drop = FALSE]
    if (!is.null(pk) && nrow(pk)) {
      p <- plotly::add_markers(p, x = pk[[tr[1]]], y = pk[[tr[2]]], z = pk$fit,
                               marker = list(size = 6, color = "#f1c40f", symbol = "diamond", line = list(color = "black", width = 1)),
                               name = "Interior peak", hoverinfo = "name")
    }
    gr <- surf$groups
    if (isTRUE(input$show_groups) && !is.null(gr) && nrow(gr) && any(!is.na(g$.fit))) {
      m1 <- gr[[paste0("mean_", tr[1])]]; m2 <- gr[[paste0("mean_", tr[2])]]
      # one label per position, so groups standardised within themselves read "a, b"
      key <- paste(round(m1, 6), round(m2, 6)); first <- !duplicated(key)
      lab <- vapply(key[first], function(k) paste(gr$group[key == k], collapse = ", "), character(1))
      m1 <- m1[first]; m2 <- m2[first]
      # each mean sits just above the nearest kept cell, so the surface does not hide it
      kept <- g[!is.na(g$.fit), , drop = FALSE]
      near <- vapply(seq_along(m1), function(i) which.min((kept[[tr[1]]] - m1[i])^2 + (kept[[tr[2]]] - m2[i])^2), integer(1))
      lift <- 0.02 * diff(range(kept$.fit))
      p <- plotly::add_trace(p, type = "scatter3d", mode = "markers+text", x = m1, y = m2, z = kept$.fit[near] + lift,
                             text = unname(lab), textposition = "top center",
                             marker = list(size = 5, color = "white", symbol = "circle", line = list(color = "black", width = 1)),
                             name = "Group mean", hoverinfo = "text")
      hp <- gr[!is.na(gr[[paste0("peak_", tr[1])]]), , drop = FALSE]
      if (nrow(hp)) {
        is_peak <- hp$peak_interior %in% TRUE
        p <- plotly::add_markers(p, x = hp[[paste0("peak_", tr[1])]], y = hp[[paste0("peak_", tr[2])]], z = hp$peak_fit,
                                 marker = list(size = 4, color = "black", symbol = ifelse(is_peak, "square", "square-open")),
                                 text = paste(hp$group, ifelse(is_peak, "peak", "high point")),
                                 name = "Group high point", hoverinfo = "text")
      }
    }
    plotly::layout(p,
      title = list(text = "Fitness surface (drag to rotate)", font = list(size = 14)),
      showlegend = FALSE, margin = list(l = 0, r = 0, b = 0, t = 40),
      scene = list(
        xaxis = list(title = sd_lab(tr[1])), yaxis = list(title = sd_lab(tr[2])),
        zaxis = list(title = z_name),
        camera = list(eye = list(x = 1.4, y = -1.4, z = 1.3))))
  })

  opt_info <- reactive({
    sf <- surfaces(); opt <- sf$land$optimum; tr <- sf$traits; g <- sf$land$grid
    dx <- opt[[tr[1]]]; dy <- opt[[tr[2]]]
    on_edge <- vapply(tr, function(t) {
      v <- g[[t]]; opt[[t]] <= min(v) + 1e-9 || opt[[t]] >= max(v) - 1e-9
    }, logical(1))
    where <- sprintf("%s = %+.2f SD and %s = %+.2f SD from the current population mean", tr[1], dx, tr[2], dy)
    sup <- sf$land$support
    # the convolution draws no population, it weighs the whole distribution
    beyond <- if (is.null(sup)) "" else sprintf(" %.0f%% of the population%s there lies outside the data.", 100 * sup$at_optimum,
                                                if (identical(sf$land$method, "convolve")) "" else " simulated")
    text <- if (any(on_edge)) {
      sprintf("No interior optimum: mean fitness keeps rising towards %s (highest %.3f).%s", where, opt$.mean_fit, beyond)
    } else {
      sprintf("Optimum at %s, %.2f SD from the current mean (mean fitness %.3f).%s", where, sqrt(dx^2 + dy^2), opt$.mean_fit, beyond)
    }
    # amber rather than green when the optimum is one to be wary of
    list(text = text, caveat = any(on_edge) || (!is.null(sup) && sup$at_optimum > 0.5))
  })
  output$opt_txt <- renderText(opt_info()$text)
  output$opt_box <- renderUI(div(class = if (opt_info()$caveat) "interp caveat" else "interp", textOutput("opt_txt")))
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
  write_groups <- function(f) {
    g <- setup()$grouped; req(!is.null(g))
    utils::write.csv(as.data.frame(g), f, row.names = FALSE)
  }
  output$dl_groups <- downloadHandler(
    filename = function() paste0(dl_name("selection_gradients_by", setup()$group), ".csv"),
    content = write_groups)
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
      theme(strip.text.y = element_text(angle = 0), strip.background = element_rect(fill = "#E6EBE7"))
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
           spline_k = band_k() %||% 10,
           surf_traits = if (length(s$traits) >= 2 && length(tr) == 2 && !anyNA(tr) && tr[1] != tr[2]) tr else character(),
           n_boot = if (!is.null(b)) attr(b, "n_drawn") else clamp_num(input$n_boot, 2, 5000, 500),
           uncertainty = if (s$surf_method == "tps" || is.null(input$surf_unc)) "none" else input$surf_unc,
           canonical = isTRUE(input$canonical), uni_land = isTRUE(input$uni_land))
  })
  output$r_code <- renderText(paste(code_lines(), collapse = "\n"))
  output$dl_code <- downloadHandler(
    filename = function() paste0(dl_name("selection_analysis"), ".R"),
    content = function(f) writeLines(code_lines(), f))

  settings_text <- reactive({
    s <- setup(); b <- boot_val()
    # mgcv turns GCV.Cp into UBRE for survival and Poisson counts and into REML
    # for a negative binomial; the fitness function takes the app's fitness
    # type, the surface finds the type in the data
    crit <- function(sm, type) if (!identical(sm, "GCV.Cp")) sm else if (type == "binary") "UBRE" else
      if (type == "count") switch(s$count_family, poisson = "UBRE", nb = "REML", "GCV") else "GCV"
    surf_type <- suppressWarnings(detect_family(s$d[[s$fit]])$type)
    surf_crit <- crit(s$surf_sm, surf_type)
    src <- DATA_SOURCES[[s$dataset]]
    paste(
      sprintf("File: %s%s", s$source, if (is.null(src)) "" else
        paste0("\n", paste(c(strwrap(src$cite, 90, initial = "Data source: ", prefix = "  "),
                             if (length(src$links)) sprintf("  %s: %s", names(src$links), src$links)),
                           collapse = "\n"))),
      sprintf("Fitness: %s (%s, %s)", s$fit, s$ftype, s$ftype_how),
      sprintf("Traits: %s", paste(s$traits, collapse = ", ")),
      if (is.null(s$group)) "Group: none" else
        sprintf("Group: %s (%s%s)", s$group,
                if (s$within_group) "standardised and relativised within group, group as a fixed effect" else
                  "one population, standardised together; group marks means and peaks on the surface",
                if (s$per_group) "; selection also estimated per group" else ""),
      sprintf("Individuals with complete data: %d", nrow(s$d)),
      sprintf("Fitness function: %s basis, %s smoothing, k = %d%s", s$spline_bs,
              crit(s$spline_sm, s$ftype), band_k(),
              if (identical(s$ftype, "count")) paste0(", ", s$count_family, " family") else ""),
      sprintf("Fitness surface: %s, %d x %d grid, %s",
              if (s$surf_method == "tps") paste0("thin-plate spline", if (s$clamp) ", held within the fitness range" else "") else
                sprintf("GAM with %s basis, %s smoothing, k %s%s", s$surf_bs, surf_crit,
                        if (is.null(s$surf_k)) "from the data" else paste("=", s$surf_k),
                        if (surf_type == "count") paste0(", ", s$count_family, " family") else ""),
              s$surf_grid, s$surf_grid, paste0(if (s$surf_full) "drawn over the full grid" else "blank outside the data",
                if (!is.null(s$surf_far)) sprintf(", cells farther than %s of the range from any individual blank", s$surf_far) else "",
                if (!is.null(s$group)) paste0(", group means and peaks marked", if (s$within_group) ", group as a fixed effect" else ", one surface for all groups") else "")),
      sprintf("Adaptive landscape: %d x %d grid, %s", s$grid_n, s$grid_n,
              switch(s$land_method,
                     convolve = "Gaussian convolution",
                     lattice = sprintf("%d simulated individuals per point, fitness interpolated from a lattice", s$sim_n),
                     sprintf("%d simulated individuals per point, fitness predicted for each", s$sim_n))),
      if (is.null(b)) "Bootstrap: not run" else sprintf("Bootstrap: %d resamples", attr(b, "n_boot")),
      sprintf("Random seed: %d", s$seed),
      sprintf("lande %s, R %s.%s", as.character(utils::packageVersion("lande")), R.version$major, R.version$minor),
      sep = "\n")
  })
  output$settings <- renderText(settings_text())

  # ---- Everything in one zip ----
  # the data, the script and every table and figure of the last run in one
  # folder, so the script runs as it is once unzipped
  output$dl_all <- downloadHandler(
    filename = function() paste0(dl_name("selection_analysis"), ".zip"),
    contentType = "application/zip",
    content = function(f) {
      s <- setup()
      dir <- file.path(tempfile("zip"), dl_name("selection_analysis"))
      dir.create(dir, recursive = TRUE)
      on.exit(unlink(dirname(dir), recursive = TRUE), add = TRUE)
      made <- character(); left_out <- character()
      # each file on its own, so one that cannot be made leaves out only itself
      add <- function(name, what, write) {
        ok <- tryCatch({ write(file.path(dir, name)); file.exists(file.path(dir, name)) }, error = function(e) FALSE)
        if (ok) made <<- c(made, sprintf("%s: %s", name, what)) else left_out <<- c(left_out, what)
      }
      figure <- function(plot, width, height) function(p) ggsave(p, plot(), width = width, height = height, dpi = 200)
      uploaded <- !s$dataset %in% names(DATA_CODE)
      withProgress(message = "Preparing the zip", value = 0, {
        if (uploaded) add(s$file_name, "the data, as uploaded", function(p) stopifnot(file.copy(s$file_path, p)))
        else add(paste0(file_stem(s$source), ".csv"), "the data, as bundled with lande",
                 function(p) utils::write.csv(get(file_stem(s$source)), p, row.names = FALSE))
        add(paste0(dl_name("selection_analysis"), ".R"), "the R code for these results",
            function(p) writeLines(code_lines(), p))
        add(paste0(dl_name("selection_gradients"), ".csv"),
            if (is.null(boot_val())) "selection differentials and gradients" else "selection differentials and gradients, with bootstrap intervals",
            write_report)
        add(paste0(dl_name("selection_gradients"), ".png"), "gradient plot", figure(grad_plot_obj, 7, 4.5))
        add(paste0(dl_name("assumption_checks"), ".csv"), "assumption checks", write_assump)
        if (!is.null(canon())) add(paste0(dl_name("canonical_axes"), ".csv"), "canonical axes of gamma", function(p) {
          ca <- canon()
          utils::write.csv(cbind(as.data.frame(ca$axes), as.data.frame(t(ca$M), check.names = FALSE)), p, row.names = FALSE)
        })
        incProgress(0.2, detail = "fitness function")
        add(paste0(dl_name("fitness_function", uni_choice()), ".png"), "fitness function", figure(uni_plot_obj, 7, 5))
        if (isTRUE(input$uni_land)) add(paste0(dl_name("adaptive_landscape", uni_choice()), ".png"),
                                        "adaptive landscape for one trait", figure(uni_land_plot_obj, 7, 5))
        if (length(s$traits) >= 2) {
          incProgress(0.3, detail = "fitness surface and landscape")
          tr <- tryCatch(surf_traits(), error = function(e) NULL)
          add(paste0(dl_name("fitness_surface", tr), ".png"), "fitness surface", figure(surf_plot_obj, 7.5, 6))
          add(paste0(dl_name("adaptive_landscape", tr), ".png"), "adaptive landscape", figure(land_plot_obj, 7.5, 6))
        }
        if (!is.null(s$grouped)) {
          add(paste0(dl_name("selection_gradients_by", s$group), ".csv"), "gradients by group", write_groups)
          add(paste0(dl_name("selection_gradients_by", s$group), ".png"), "gradient plot by group", figure(group_plot_obj, 7, 4.5))
        }
        incProgress(0.4, detail = "zipping")
        writeLines(c(
          "lande selection analysis",
          format(Sys.Date(), "%d %B %Y"), "",
          "Files", paste0("  ", c(made, "README.txt: this file")), "",
          if (length(left_out)) c("Left out, as they could not be made with these settings", paste0("  ", left_out), ""),
          "Settings", paste0("  ", strsplit(settings_text(), "\n", fixed = TRUE)[[1]]), "",
          "To repeat the analysis",
          "  Install lande, set R's working directory to this folder (in RStudio: Session >",
          "  Set Working Directory > To Source File Location) and run the R script.",
          if (uploaded) "  It reads the data file beside it." else "  It loads the data from lande; the CSV is a copy.",
          "  The figures here follow the display options set in the app, so the plots the",
          "  script draws can look different."
        ), file.path(dir, "README.txt"))
        zip_folder(f, dir)
      })
    })

  # the maintenance notice (see read_notice): the file is looked at every 30
  # seconds, which costs one stat when it is absent, and read only when it
  # changes. Sent ahead of the other observers, so a page gets it at once.
  if (nzchar(NOTICE_FILE)) {
    notice <- reactivePoll(30000, session,
      checkFunc = function() paste(unlist(file.info(NOTICE_FILE, extra_cols = FALSE)[c("size", "mtime", "ctime")]), collapse = " "),
      valueFunc = function() read_notice(NOTICE_FILE))
    observe(session$sendCustomMessage("lande_notice", notice()), priority = 100)
  }
}

shinyApp(ui, server)
