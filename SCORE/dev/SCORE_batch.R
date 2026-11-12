# =============================================================================
# SCORE_batch.R
# Runs the SCORE workflow (Tabs 1-4) outside the Shiny app, for a batch of
# pilots, and saves every run exactly like save_run() did (parameters, best
# schedule, all candidates, tables, plots, index).
#
# What is reproduced from SCORE.R (same functions from functions.R, same rules):
#   Tab 1  import: timestamp from Begin.Path (first YYYYMMDD_HHMMSS of the file
#          name), date = date of the file start, minute = clock minute of the
#          detection, recorder "global" (or parsed from Begin.Path)
#   Tab 2  duty-cycle schedules (create_sun_duty_schedule), bootstrap over days,
#          rarefaction curve and missed-species heatmap
#   Tab 3  time-window grid W1 / W2 / combined (top N per duration), richness
#          normalised to the best design
#   Tab 4  multi-optimum: every duty cycle x every W1, W2 and top-N combined
#          window, blocks kept with block_in_window() (corrected version),
#          mean daily richness, Pareto front
#
# Speed: Tab 4 evaluates all candidates with matrix operations instead of one
# compute_richness() call per candidate. Each run then RE-COMPUTES a random
# sample of candidates (and the selected one) with the app code path
# (schedule rows + block_in_window() + compute_richness()) and stops if any
# value differs. So the speed-up cannot change the results silently.
#
# Three chapters, which can be run separately (switches in section 0a):
#   CHAPTER A  Temporal pilots: select on each pilot (start date x length) and
#              save the run (then validate with validate_pilots.R)
#   CHAPTER B  Spatial validation: on the same days, select on a subset of the
#              recorders (training) and evaluate the selected schedule on the
#              other recorders (test), relative to continuous 24-h recording
#   CHAPTER C  SCORE presets: runs Tabs 2-4 on chosen periods and saves a preset
#              (.zip with the .json and the figures), identical to the one of the
#              app tab "6. Share results" (same code: preset_functions.R)
#
# Not included: Tab 5 (target species), which is not used for the pilots.
# =============================================================================

## ---- 0. SETTINGS (edit here) ------------------------------------------------

## 0a. What to run -------------------------------------------------------------
run_temporal <- FALSE     # CHAPTER A: temporal pilots
run_spatial  <- FALSE     # CHAPTER B: spatial validation
run_presets  <- TRUE      # CHAPTER C: SCORE presets (.zip, as in the app)

## 0b. CHAPTER A: pilot start dates x pilot lengths (12 pilots as planned) -----
pilot_starts  <- c("2023-05-12", "2023-06-01", "2023-07-01")
pilot_lengths <- c(1, 3, 7, 14)
# To run a custom list instead, give a table, e.g.:
# pilots_custom <- data.frame(start = c("2023-06-01", "2023-08-01"), length = c(7, 3))
pilots_custom <- NULL

## 0c. CHAPTER B: spatial validation --------------------------------------------
spatial_start   <- "2023-06-01"     # same days for training and test recorders
spatial_end     <- "2023-06-07"
spatial_splits  <-  c("group", "random")
#   "group"  = leave-one-group-out: select on one group (e.g. site), test on the others
#   "random" = random splits: n_train recorders for selection, the others for testing
#   "loo"    = leave-one-out: select on all recorders but one, test on that one
unit_group      <- function(unit) sub("[0-9]+$", "", unit)   # recorder -> group: GONDO1 -> GONDO
n_random_splits <- 7
n_train         <- 10               # recorders used for selection in each random split
# Recorders ("units") are identified by the file-name prefix (GONDO1_20230601_050000 -> GONDO1).
# If the recorder is a folder of Begin.Path instead, give its V index (as in the app's path parser):
unit_path_col   <- NULL

## 0c'. CHAPTER C: presets -------------------------------------------------------
# Fields of the app form, common to all presets (each preset below can override any of them)
preset_common <- list(
  name = "Gabriel Marcacci", affiliation = "Swiss Ornithological Institute", orcid = "", contact = "gabriel.marcacci@vogelwarte.ch; amandine.serrurier@vogelwarte.ch",
  project_description = "17 in Burkina Faso 24h for 3 days",
  objective      = "community",            # or "target species", or your own text
  target_taxa    = "birds",
  target_species = NULL,                   # e.g. c("Rock Ptarmigan", "Alpine Accentor")
  country        = "Burkina Faso",          # full name, as in the app list (COUNTRIES)
  region         = "",
  habitat_codes  = c("E", "G"),                 # EUNIS level 1 codes, e.g. c("E", "H")
  habitat_desc   = "Sudanian savanna",
  elevation_mean = 200,                     # m
  share_location = TRUE,  round_location = TRUE,   # coordinates of 0e, rounded to ~10 km
  timing         = "December",                   # e.g. c("breeding", "June")
  recorder_model = "SMMicro", sample_rate = "", height = "",
  raw_recording  = "continuous 24 h",
  classifier     = "BirdNET-Analyzer", classifier_version = "2.4",
  notes          = "no confidence threshold")
# One element per preset: period + fields that differ from preset_common.
# Optional: s_override = list(...) changes settings of section 0e-0h for this preset only
#           (e.g. list(period_step = 1) to also test 1/2 and 1/3 min; ~3x slower than step 3)
# Optional: validation = list(...) is copied into the preset (shown on the GitHub page)

presets_to_make <- list(
  list(start = "2025-12-31", end = "2026-01-02",          # tes 3 jours de décembre 2025
       title = "Sudanian savanna bird community, Burkina Faso: December")
)
# presets_to_make <- list(
#   list(start = "2023-06-01", end = "2026-06-07",
#        title = "Burkina Faso community, 3 days, December",
#        timing = c("December"),
#        notes = "Reference preset from Serrurier et al. (JAV, in revision); same settings as in the paper.",
#        validation = list(      # results of validate_pilots.R and CHAPTER B for this pilot (paper settings)
#          temporal = list(description = "Schedule selected on this pilot (90% target) applied to other days of the same season, same recorders",
#                          in_sample_pct_of_24h = 89.7, later_days_pct_of_24h = 78.2, later_days_n = 99,
#                          all_other_days_pct_of_24h = 79.5),
#          spatial = list(description = "Schedules selected on subsets of these recorders (same 7 days) and applied to the other recorders",
#                         between_sites = list(gondo_to_grengiols = 86.7, grengiols_to_gondo = 86.5),
#                         random_10_to_4 = list(mean = 75.6, min = 69.6, max = 82.3, n_runs = 7),
#                         leave_one_recorder_out = list(mean = 59.2, min = 46.5, max = 71.6, n_runs = 14)),
#          guidance = "Transfers well within the same season phase and between these two Alpine sites; not to single recorders or to late summer (69-83% from late July).")),
#   list(start = "2023-07-01", end = "2023-07-07",
#        title = "Alpine bird community, Valais (CH): late breeding season",
#        timing = c("breeding", "July")),
#   list(start = "2023-08-25", end = "2023-08-31",
#        title = "Alpine bird community, Valais (CH): post-breeding season",
#        timing = c("post-breeding", "August"))
# )
preset_content <- list(targets = c(75, 80, 85, 90, 95), include_pareto = TRUE,
                       include_species = TRUE, species_min_conf = 0)
SCORE_VERSION <- "1.1.0"

## 0d. Files ------------------------------------------------------------------
score_dir  <-  "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/VoWa/shiny_App_raref/bioacoustic-tools/SCORE"   # folder with functions.R
data_files <- "/Users/amandine.serrurier/Documents/traning_data_complete_bn_2024_GG/BurkinaGM_2025-12/all_data_filtered_BF_1225.txt"
runs_dir         <- "/Users/amandine.serrurier/Documents/SCORE_runs/batch_runs"   # CHAPTER A, one sub-folder per run
spatial_runs_dir <- "/Users/amandine.serrurier/Documents/SCORE_runs/spatial_runs" # CHAPTER B
presets_out_dir  <- "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/VoWa/shiny_App_raref/bioacoustic-tools/SCORE/presets_library"
library_presets_dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/VoWa/shiny_App_raref/bioacoustic-tools/SCORE/presets_library/presets"  # NULL = ne rien copier# CHAPTER C: one .zip per preset
skip_existing_runs <- TRUE      # TRUE = a selection (or preset .zip) already saved is not recomputed

## 0e. Tab 1: import, filters, coordinates (same meaning as in the app) -------
s <- list(
  timezone          = "Africa/Ouagadougou",
  min_conf          = 0,
  species_file      = NULL,   
  exclude_recorders = NULL,
  site_col          = NULL,
  rec_col           = NULL,
  coord_mode        = "single",
  single_lat        = 12.25,    
  single_lon        = -0.63,    


## 0f. Tab 2: duty cycle ----------------------------------------------------
  record_24h        = TRUE,
  start_event       = "dawn",   end_event   = "dawn",          # used only if record_24h = FALSE
  length_morning    = 0,        extend_before = 1, extend_after = 0,
  use_window2       = FALSE,
  start_event2      = "dawn",   end_event2  = "dawn",
  length_morning2   = 0,        extend_before2 = 1, extend_after2 = 0,
  period_min        = 1,        period_max  = 60,  period_step = 3,
  duty_duration     = 1,
  B_temporal        = 200,      # bootstrap replicates (Tab 2 plots only)
  skip_bootstrap_1day = TRUE,   # 1-day pilot: nothing to resample (see note in section 5)
  show_summary      = TRUE,

## 0g. Tab 3: time windows --------------------------------------------------
  grid_event        = "sunrise",
  grid_start_min    = -3,   grid_start_max  = 5,  grid_start_step  = 0.5,
  grid_dur_min      = 0.5,  grid_dur_max    = 8,  grid_dur_step    = 0.5,
  use_window2_grid  = TRUE,
  grid_event2       = "sunset",
  grid_start_min2   = -5,   grid_start_max2 = 2,  grid_start_step2 = 0.5,
  grid_dur_min2     = 0.5,  grid_dur_max2   = 8,  grid_dur_step2   = 0.5,
  top_n_combined    = 2,
  show_frontier     = TRUE,

## 0h. Tab 4 and saving -----------------------------------------------------
  threshold         = 90,       # best schedule = least effort with richness_pct >= threshold
  near_optimal_pct  = 5,        # also save schedules within +5 % of the best effort
  n_verify          = 25,       # candidates re-computed with the app code path in each run
  seed              = 1,        # random seed (bootstrap, verification sample)
  save_plotly       = TRUE      # also save interactive .html versions of the plots
)
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(lubridate); library(suncalc); library(dplyr)
  library(stringr); library(ggplot2)
})
has_scico  <- requireNamespace("scico", quietly = TRUE)
has_plotly <- requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)
# need_cols() in functions.R uses shiny's validate()/need()
if (requireNamespace("shiny", quietly = TRUE)) { validate <- shiny::validate; need <- shiny::need } else {
  need <- function(expr, message) if (isTRUE(expr)) NULL else message
  validate <- function(...) { m <- unlist(list(...)); if (length(m)) stop(m, call. = FALSE); invisible(TRUE) }
}
source(file.path(score_dir, "functions.R"))
if (run_presets) source(file.path(score_dir, "preset_functions.R"))   # same code as the app (Tab 6)
if (!isTRUE(all.equal(block_in_window(1:3, 1:3, c(0, 2, 2), c(1, 3, 3)), c(TRUE, TRUE, TRUE))))
  stop("functions.R still contains the ORIGINAL block_in_window(). Use the corrected version.")

battery_specs <- data.frame(                       # same table as in SCORE.R
  recorder     = c("SMMicro 2", "SMMini 2", "SMMini 2", "SM4", "SM5", "AudioMoth -v1.2.0"),
  nb_battery   = c(4, 8, 6, 4, 11, 3),
  battery_type = c("AA", "AA", "lithium", "4D alkaline", "lithium", "AA"),
  duration_h   = c(280, 530, 1330, 310, 1710, 189))
battery_specs$code <- paste0(gsub(" ", "", battery_specs$recorder), "_",
                             battery_specs$nb_battery, substr(battery_specs$battery_type, 1, 2))

msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
t_start <- Sys.time()

## ---- 1. Tab 1: import (once for all pilots) ---------------------------------
msg("Reading data ...")
raw <- rbindlist(lapply(data_files, fread, showProgress = FALSE), fill = TRUE)
names(raw) <- make.names(names(raw))
need_cols(raw, c("Begin.Path", "Common.Name", "Confidence", "Begin.Time..s."), "BirdNET table")
raw[, Begin.Path := normalize_begin_path(Begin.Path)]

# Same parsing as dt_parsed() in SCORE.R, written with integer arithmetic for speed
raw[, filename := tools::file_path_sans_ext(basename(Begin.Path))]
uf  <- unique(raw$filename)
hit <- regexpr("\\d{8}_\\d{6}", uf)
uts <- ifelse(hit > 0, substr(uf, hit, hit + 14L), NA_character_)
raw[, timestamp := uts[match(filename, uf)]]
raw <- raw[!is.na(timestamp)]
u <- unique(raw[, .(timestamp)])
u[, `:=`(fdate = as.Date(substr(timestamp, 1, 8), "%Y%m%d"),
         sod   = as.integer(substr(timestamp, 10, 11)) * 3600L +
                 as.integer(substr(timestamp, 12, 13)) * 60L + as.integer(substr(timestamp, 14, 15)))]
raw <- u[raw, on = "timestamp"]
raw[, start := as.numeric(Begin.Time..s.)]
raw[, minute_of_day := as.integer(((sod + floor(start)) %/% 60) %% 1440)]

# Check the fast parsing against the app code (parse_date_time + hour/minute) on a sample
set.seed(s$seed)
chk <- raw[sample.int(.N, min(.N, 50000))]
chk[, ts_adj := parse_date_time(timestamp, orders = "Ymd_HMS", tz = s$timezone) + seconds(start)]
n_bad <- chk[minute_of_day != hour(ts_adj) * 60 + minute(ts_adj), .N]
if (n_bad > 0) stop(n_bad, " detections get a different minute than in the app (DST change in the data?).")
rm(chk)

dtA <- raw[, .(species = Common.Name, conf = Confidence, date = fdate, minute_of_day,
               filename, Begin.Path)]
rm(raw, u); invisible(gc())
dtA[, recorder_file := sub("_\\d{8}_\\d{6}.*$", "", filename)]
if (length(s$exclude_recorders)) dtA <- dtA[!recorder_file %in% s$exclude_recorders]
# Recorder identity kept as "unit" (used only to split recorders in CHAPTER B; the analysis
# itself pools all recorders per day, as in the app's "global" mode)
dtA[, unit := if (is.null(unit_path_col)) recorder_file else split_path_matrix(Begin.Path)[, unit_path_col]]

# Path parsing (as "Apply path parsing" in the app)
if (!is.null(s$site_col) || !is.null(s$rec_col)) {
  mat <- split_path_matrix(dtA$Begin.Path)
  dtA[, site     := if (!is.null(s$site_col)) mat[, s$site_col] else "all"]
  dtA[, recorder := if (!is.null(s$rec_col))  mat[, s$rec_col]  else "global"]
  rm(mat)
} else {
  dtA[, `:=`(site = "all", recorder = "global")]
}
# Filters (as "Apply filters": min confidence, species)
dtA <- dtA[conf >= s$min_conf]
if (!is.null(s$species_file)) {
  ext <- tolower(tools::file_ext(s$species_file))
  sp_keep <- if (ext == "txt") {
    x <- trimws(sub("\r$", "", readLines(s$species_file, encoding = "UTF-8", warn = FALSE)))
    unique(sub("^[^_]*_", "", x[x != ""]))
  } else if (ext %in% c("xlsx", "xls")) readxl::read_excel(s$species_file)$species else read.csv(s$species_file)$species
  dtA <- dtA[species %in% trimws(as.character(sp_keep))]
}
dtA[, c("filename", "Begin.Path", "recorder_file") := NULL]
setkey(dtA, date)
species_levels <- sort(unique(dtA$species))
msg(format(nrow(dtA), big.mark = ","), " detections, ", length(species_levels), " species, ",
    uniqueN(dtA$date), " days, recorders: ", paste(unique(dtA$recorder), collapse = ", "))

# Coordinates (as in the app)
coords <- if (s$coord_mode == "multiple") read_coord_file(s$coord_file, s$coord_file) else NULL

solar_for <- function(rec_dates) {                 # rec_dates: data.table(recorder, date)
  if (s$coord_mode == "single") {
    d <- sort(unique(rec_dates$date))
    sun <- as.data.table(getSunlightTimes(date = d, lat = s$single_lat, lon = s$single_lon,
                                          keep = c("dawn", "sunrise", "sunset", "dusk"), tz = s$timezone))
    sun[, date := as.Date(date)]
    out <- merge(unique(rec_dates[, .(recorder, date)]), sun[, .(date, dawn, sunrise, sunset, dusk)], by = "date")
  } else {
    comb <- merge(unique(rec_dates[, .(recorder, date)]), coords, by = "recorder", all.x = TRUE)
    out <- rbindlist(lapply(seq_len(nrow(comb)), function(i) {
      x <- getSunlightTimes(date = as.Date(comb$date[i]), lat = comb$latitude[i], lon = comb$longitude[i],
                            keep = c("dawn", "sunrise", "sunset", "dusk"), tz = s$timezone)
      data.table(recorder = comb$recorder[i], date = as.Date(comb$date[i]),
                 dawn = x$dawn, sunrise = x$sunrise, sunset = x$sunset, dusk = x$dusk)
    }))
  }
  setcolorder(out, c("recorder", "date")); out[]
}
ref_minute <- function(sol, event) {
  if (event %in% c("dawn", "sunrise", "sunset", "dusk"))
    as.integer(hour(sol[[event]]) * 60 + minute(sol[[event]]))
  else rep(as.integer(sub(":.*", "", event)) * 60L + as.integer(sub(".*:", "", event)), nrow(sol))
}

## ---- 2. Tab 2: duty-cycle schedules and bootstrap -------------------------
run_tab2 <- function(dt, sol_all, bootstrap = TRUE) {
  periods    <- seq(as.numeric(s$period_min), as.numeric(s$period_max), by = as.numeric(s$period_step))
  param_grid <- CJ(period = periods, duty_duration = s$duty_duration)
  recs <- unique(dt$recorder)
  schedules <- lapply(seq_len(nrow(param_grid)), function(i) {
    p <- param_grid[i]
    sched <- rbindlist(lapply(recs, function(r) {
      sun_r <- sol_all[recorder == r]
      s1 <- create_sun_duty_schedule(sun_table = sun_r, start_event = s$start_event, end_event = s$end_event,
                                     length_morning = s$length_morning, extend_before_hours = s$extend_before,
                                     extend_after_hours = s$extend_after, period = p$period,
                                     duty_duration = p$duty_duration, nb_duty = 1, random_start = FALSE,
                                     record_24h = s$record_24h)
      if (isTRUE(s$use_window2)) {
        s2 <- create_sun_duty_schedule(sun_table = sun_r, start_event = s$start_event2, end_event = s$end_event2,
                                       length_morning = s$length_morning2, extend_before_hours = s$extend_before2,
                                       extend_after_hours = s$extend_after2, period = p$period,
                                       duty_duration = p$duty_duration, nb_duty = 1, random_start = FALSE)
        s1 <- rbindlist(list(s1, s2))
      }
      s1[, recorder := r]; s1
    }))
    setattr(sched, "period", p$period); setattr(sched, "duty_duration", p$duty_duration)
    sched
  })

  if (!bootstrap) return(list(schedules = schedules, summary = NULL, missed = NULL, B = 0L))
  all_dates <- unique(dt$date); n_days <- length(all_dates)
  filtered_dates_dt <- data.table(date = all_dates)
  schedule_periods  <- sapply(schedules, attr, "period")
  schedule_efforts  <- sapply(schedules, function(x) {
    x2 <- merge(as.data.table(x), filtered_dates_dt, by = "date"); sum(x2$end_min - x2$start_min + 1) })

  # Daily species detected under each schedule (= what compute_richness() returns per day)
  Dd <- unique(dt[, .(recorder, date, minute = as.numeric(minute_of_day), species)])
  daily <- lapply(schedules, function(x) {
    x <- as.data.table(x)[date %in% all_dates][, `:=`(start_min = as.numeric(start_min), end_min = as.numeric(end_min))]
    h <- x[Dd, on = .(recorder, date, start_min <= minute, end_min >= minute), nomatch = 0L,
           .(date = i.date, species = i.species)]
    unique(h)
  })

  B <- if (n_days == 1 && isTRUE(s$skip_bootstrap_1day)) 0L else as.integer(s$B_temporal)
  if (B == 0L) return(list(schedules = schedules, summary = NULL, missed = NULL, B = 0L))

  all_results <- vector("list", B)
  for (b in seq_len(B)) {
    # NOTE: the app uses sample(all_dates, ...), which goes wrong when there is only one
    # date (sample() then draws from 1:x); indexing avoids that and is otherwise identical.
    sampled_days <- all_dates[sample.int(n_days, n_days, replace = TRUE)]
    one <- lapply(seq_along(schedules), function(i) {
      d <- daily[[i]]
      r_day <- d[, .(richness = uniqueN(species)), by = date]
      r_m <- merge(data.table(date = sampled_days), r_day, by = "date", all.x = TRUE)
      r_m[is.na(richness), richness := 0]
      data.table(period = schedule_periods[i], duty_duration = attr(schedules[[i]], "duty_duration"),
                 effort = schedule_efforts[i], richness = mean(r_m$richness), bootstrap = b,
                 species_detected = list(unique(d[date %in% sampled_days, species])))
    })
    all_results[[b]] <- rbindlist(one)
  }
  summary_richness <- rbindlist(all_results)
  summary_richness[, duty_cycle_label := paste0(duty_duration, "/", period)]
  summary_richness[, richness_pct := 100 * richness / max(richness), by = bootstrap]
  summary_richness[, n_species_detected := lengths(species_detected)]

  ## Missed-species analysis (copied from SCORE.R)
  min_period_val  <- min(summary_richness$period)
  reference_table <- summary_richness[period == min_period_val,
                                      list(reference_species = list(unique(unlist(species_detected)))), by = "bootstrap"]
  summary_richness <- merge(summary_richness, reference_table, by = "bootstrap", all.x = TRUE)
  missed_species_list <- summary_richness[, list(species_missed = mapply(
    function(ref, detected) setdiff(ref, unlist(detected)), reference_species, species_detected, SIMPLIFY = FALSE)),
    by = list(bootstrap, period)]
  missed_species_long <- missed_species_list[, list(species_missed = unlist(species_missed)), by = list(bootstrap, period)]
  missed_summary <- missed_species_long[, .N, by = list(species_missed, period)]
  species_order <- missed_summary[, list(total_missed = sum(N)), by = "species_missed"][order(total_missed), species_missed]
  all_periods <- sort(unique(summary_richness$period))
  missed_summary_full <- CJ(species_missed = species_order, period = all_periods)
  missed_summary_full <- merge(missed_summary_full, missed_summary, by = c("species_missed", "period"), all.x = TRUE)
  missed_summary_full[is.na(N), N := 0]
  missed_summary_full[, species_missed := factor(species_missed, levels = species_order)]
  list(schedules = schedules, summary = summary_richness, missed = missed_summary_full, B = B)
}

## ---- 3. Tab 3: time-window grid (code from SCORE.R) -------------------------
run_tab3 <- function(dt_in, sol_all) {
  dt <- unique(dt_in[, .(recorder, date, minute = minute_of_day, species)])   # same species sets
  setkey(dt, recorder, date)
  comb <- unique(dt[, .(recorder, date)])
  comb_solar <- sol_all[, .(recorder, date, dawn, sunrise, sunset, dusk)][comb, on = c("recorder", "date")]
  start_offsets_min <- seq(s$grid_start_min * 60, s$grid_start_max * 60, by = s$grid_start_step * 60)
  durations_min     <- seq(s$grid_dur_min * 60,   s$grid_dur_max * 60,   by = s$grid_dur_step * 60)
  comb_solar[, ref_min := ref_minute(comb_solar, s$grid_event)]
  if (isTRUE(s$use_window2_grid)) comb_solar[, ref_min2 := ref_minute(comb_solar, s$grid_event2)]

  species_in_window <- function(rec_i, d_i, s0_i, s1_i) {
    if (s0_i >= 0 && s1_i <= 1440) {
      unique(dt[.(rec_i, d_i)][minute >= s0_i & minute <= s1_i, species])
    } else if (s1_i > 1440) {
      sp1 <- unique(dt[.(rec_i, d_i)][minute >= max(0L, s0_i) & minute <= 1440L, species])
      sp2 <- unique(dt[.(rec_i, as.Date(d_i) + 1L)][minute >= 0L & minute <= (s1_i - 1440L), species])
      unique(c(sp1, sp2))
    } else if (s0_i < 0) {
      sp0 <- unique(dt[.(rec_i, as.Date(d_i) - 1L)][minute >= (1440L + s0_i) & minute <= 1440L, species])
      sp1 <- unique(dt[.(rec_i, d_i)][minute >= 0L & minute <= min(1440L, s1_i), species])
      unique(c(sp0, sp1))
    } else character(0)
  }

  grid_all <- comb_solar[, as.data.table(expand.grid(start_offset_min = start_offsets_min,
                                                     duration = durations_min)), by = .(recorder, date, ref_min)]
  grid_all[, `:=`(start_abs_min = ref_min + start_offset_min, effort = duration)]
  grid_all[, species_w1 := mapply(species_in_window, recorder, date, start_abs_min,
                                  start_abs_min + duration, SIMPLIFY = FALSE)]
  grid_all[, richness := sapply(species_w1, length)]
  grid_all[, start_off_h := start_offset_min / 60]
  grid_all[, duration_h  := duration / 60]
  grid_daily <- grid_all[, .(richness = length(unique(unlist(species_w1))), effort = mean(effort, na.rm = TRUE),
                             species_union = list(unique(unlist(species_w1)))), by = .(date, start_off_h, duration_h)]
  grid_w1 <- grid_daily[, .(richness = mean(richness, na.rm = TRUE), effort = mean(effort, na.rm = TRUE),
                            n_species_detected = length(unique(unlist(species_union)))), by = .(start_off_h, duration_h)]
  grid_w2 <- NULL; grid_combined <- NULL

  if (isTRUE(s$use_window2_grid)) {
    start_offsets_min2 <- seq(s$grid_start_min2 * 60, s$grid_start_max2 * 60, by = s$grid_start_step2 * 60)
    durations_min2     <- seq(s$grid_dur_min2 * 60,   s$grid_dur_max2 * 60,   by = s$grid_dur_step2 * 60)
    grid_all2 <- comb_solar[, as.data.table(expand.grid(start_offset_min = start_offsets_min2,
                                                        duration = durations_min2)), by = .(recorder, date, ref_min = ref_min2)]
    grid_all2[, `:=`(start_abs_min = ref_min + start_offset_min, effort = duration)]
    grid_all2[, species_w2 := mapply(species_in_window, recorder, date, start_abs_min,
                                     start_abs_min + duration, SIMPLIFY = FALSE)]
    grid_all2[, richness := sapply(species_w2, length)]
    grid_all2[, start_off_h := start_offset_min / 60]
    grid_all2[, duration_h  := duration / 60]
    grid_daily2 <- grid_all2[, .(richness = length(unique(unlist(species_w2))), effort = mean(effort, na.rm = TRUE),
                                 species_union = list(unique(unlist(species_w2)))), by = .(date, start_off_h, duration_h)]
    grid_w2 <- grid_daily2[, .(richness = mean(richness, na.rm = TRUE), effort = mean(effort, na.rm = TRUE),
                               n_species_detected = length(unique(unlist(species_union)))), by = .(start_off_h, duration_h)]
    top_n <- as.integer(s$top_n_combined)
    best_w1_keys <- grid_all[, .(richness = length(unique(unlist(species_w1)))), by = .(start_off_h, duration_h)][
      order(duration_h, -richness)][, .SD[seq_len(min(.N, top_n))], by = duration_h]
    best_w2_keys <- grid_all2[, .(richness = length(unique(unlist(species_w2)))), by = .(start_off_h, duration_h)][
      order(duration_h, -richness)][, .SD[seq_len(min(.N, top_n))], by = duration_h]
    grid_all_sp  <- grid_all[best_w1_keys, on = c("start_off_h", "duration_h"), nomatch = 0][
      , .(recorder, date, start_off_h, duration_h, dur_w1 = duration, species_w1)]
    grid_all2_sp <- grid_all2[best_w2_keys, on = c("start_off_h", "duration_h"), nomatch = 0][
      , .(recorder, date, start_off_h2 = start_off_h, duration_h2 = duration_h, dur_w2 = duration, species_w2)]
    grid_all_comb <- merge(grid_all_sp, grid_all2_sp, by = c("recorder", "date"), allow.cartesian = TRUE)
    grid_all_comb[, effort := dur_w1 + dur_w2]
    grid_comb_daily <- grid_all_comb[, .(
      richness = length(unique(c(unlist(species_w1), unlist(species_w2)))),
      richness_w1 = length(unique(unlist(species_w1))), richness_w2 = length(unique(unlist(species_w2))),
      effort = mean(effort, na.rm = TRUE),
      species_union = list(unique(c(unlist(species_w1), unlist(species_w2))))),
      by = .(date, start_off_h, duration_h, start_off_h2, duration_h2)]
    grid_combined <- grid_comb_daily[, .(
      richness = mean(richness, na.rm = TRUE), richness_w1 = mean(richness_w1, na.rm = TRUE),
      richness_w2 = mean(richness_w2, na.rm = TRUE), effort = mean(effort, na.rm = TRUE),
      n_species_detected = length(unique(unlist(species_union)))),
      by = .(start_off_h, duration_h, start_off_h2, duration_h2)]
    max_global <- max(c(grid_w1$richness, grid_w2$richness, grid_combined$richness), na.rm = TRUE)
    grid_w1[, richness_pct := 100 * richness / max_global]
    grid_w2[, richness_pct := 100 * richness / max_global]
    grid_combined[, richness_pct    := 100 * richness    / max_global]
    grid_combined[, richness_w1_pct := 100 * richness_w1 / max_global]
    grid_combined[, richness_w2_pct := 100 * richness_w2 / max_global]
  } else {
    grid_w1[, richness_pct := 100 * richness / max(grid_w1$richness, na.rm = TRUE)]
  }
  list(w1 = grid_w1, w2 = grid_w2, comb = grid_combined)
}

## ---- 4. Tab 4: multi-optimum --------------------------------------------------
window_sets_from <- function(grid) {
  ws <- if (!is.null(grid$w2)) list(W1 = grid$w1, W2 = grid$w2, Combined = grid$comb) else list(W1 = grid$w1)
  if (!is.null(ws$Combined)) {
    top_n <- as.integer(s$top_n_combined)
    ws$Combined <- ws$Combined[order(duration_h, -richness_pct)][, .SD[seq_len(min(.N, top_n))],
                                                                 by = .(duration_h, duration_h2)]
  }
  ws
}
cand_table <- function(window_sets) {           # same order as the app loops: set > window
  rbindlist(lapply(names(window_sets), function(nm) {
    g <- window_sets[[nm]]
    if (nm == "W1") data.table(window_type = nm, start1_hour = g$start_off_h, dur1_h = g$duration_h,
                               start2_hour = NA_real_, dur2_h = NA_real_)
    else if (nm == "W2") data.table(window_type = nm, start1_hour = NA_real_, dur1_h = NA_real_,
                                    start2_hour = g$start_off_h, dur2_h = g$duration_h)
    else data.table(window_type = nm, start1_hour = g$start_off_h, dur1_h = g$duration_h,
                    start2_hour = g$start_off_h2, dur2_h = g$duration_h2)
  }))[, w := .I][]
}

# Fast engine: for each duty schedule and each recorder-day, a species x block incidence
# matrix M and a block x candidate matrix K (block_in_window); detected = (M %*% K) > 0
run_tab4 <- function(dt, sol_all, schedules, cands) {
  pilot_dates <- sort(unique(dt$date))
  refs <- sol_all[, .(recorder, date)]
  refs[, ref_min_w1 := ref_minute(sol_all, s$grid_event)]
  refs[, ref_min_w2 := if (isTRUE(s$use_window2_grid)) ref_minute(sol_all, s$grid_event2) else ref_min_w1]
  Dd <- unique(dt[, .(recorder, date, minute = as.numeric(minute_of_day), sp = match(species, species_levels))])
  nC <- nrow(cands)
  has1 <- !is.na(cands$start1_hour); has2 <- !is.na(cands$start2_hour)
  out <- vector("list", length(schedules))
  for (p in seq_along(schedules)) {
    sched <- merge(as.data.table(schedules[[p]]), refs, by = c("recorder", "date"), all.x = TRUE)
    sched <- sched[date %in% pilot_dates]
    rich_sum <- numeric(nC); n_kept_days <- integer(nC); effort <- numeric(nC)
    sp_any <- matrix(FALSE, length(species_levels), nC)
    for (di in seq_along(pilot_dates)) {
      d <- pilot_dates[di]                            # (a for loop over Dates drops the class)
      day_det <- matrix(FALSE, length(species_levels), nC); day_kept <- rep(FALSE, nC)
      for (r in unique(sched[date == d, recorder])) {
        b <- sched[date == d & recorder == r]
        nb <- nrow(b); if (nb == 0) next
        r1 <- b$ref_min_w1[1]; r2 <- b$ref_min_w2[1]
        bi <- function(ws, we) { n <- length(ws)
          matrix(block_in_window(rep(b$start_min, n), rep(b$end_min, n),
                                 rep(ws, each = nb), rep(we, each = nb)), nb, n) }
        K <- matrix(FALSE, nb, nC)
        if (any(has1)) K[, has1] <- bi(r1 + cands$start1_hour[has1] * 60, r1 + (cands$start1_hour[has1] + cands$dur1_h[has1]) * 60)
        if (any(has2)) K[, has2] <- K[, has2] | bi(r2 + cands$start2_hour[has2] * 60, r2 + (cands$start2_hour[has2] + cands$dur2_h[has2]) * 60)
        effort <- effort + colSums(K * (b$end_min - b$start_min + 1))
        day_kept <- day_kept | colSums(K) > 0
        h <- Dd[date == d & recorder == r]
        if (nrow(h)) {
          bk <- b[, .(bstart = as.numeric(start_min), bend = as.numeric(end_min), bidx = .I)]
          inc <- unique(bk[h, on = .(bstart <= minute, bend >= minute), nomatch = 0L, .(bidx, sp = i.sp)])
          if (nrow(inc)) {
            M <- matrix(0, length(species_levels), nb); M[cbind(inc$sp, inc$bidx)] <- 1
            day_det <- day_det | (M %*% K) > 0
          }
        }
      }
      rich_sum    <- rich_sum + colSums(day_det) * day_kept
      n_kept_days <- n_kept_days + day_kept
      sp_any <- sp_any | (day_det & rep(day_kept, each = nrow(day_det)))
    }
    keep <- effort > 0                                # the app skips schedules with no block
    out[[p]] <- data.table(w = cands$w[keep], p = p, duty_period = attr(schedules[[p]], "period"),
                           effort = effort[keep], richness = rich_sum[keep] / n_kept_days[keep],
                           n_species_detected = colSums(sp_any)[keep])
  }
  res <- rbindlist(out)
  res <- cands[res, on = "w"]
  setorder(res, w, p)                                 # app order: window outer, duty cycle inner
  res[, richness_pct := 100 * richness / max(richness, na.rm = TRUE)]
  setorder(res, effort)
  # Pareto dominance, as recomputed in the app plot (gg_optimum), vectorised
  o  <- order(res$effort); e <- res$effort[o]; rr <- res$richness_pct[o]
  grp_first <- match(e, e)                           # first row with the same effort
  cm <- cummax(rr); prev_max <- c(-Inf, cm)[grp_first]   # best richness at strictly lower effort
  same_max <- ave(rr, e, FUN = max)
  dom <- logical(length(o)); dom[o] <- (prev_max >= rr) | (same_max > rr)
  res[, `:=`(pareto = !dom, dominated = dom)]
  res[, c("w", "p") := NULL]
  setcolorder(res, c("window_type", "start1_hour", "dur1_h", "start2_hour", "dur2_h", "duty_period",
                     "effort", "richness", "n_species_detected", "richness_pct", "pareto", "dominated"))
  res[]
}

# Exact app code path for one candidate (used to verify the fast engine)
app_candidate <- function(dt, sol_all, schedules, row) {
  p <- match(row$duty_period, sapply(schedules, attr, "period"))
  refs <- sol_all[, .(recorder, date)]
  refs[, ref_min_w1 := ref_minute(sol_all, s$grid_event)]
  refs[, ref_min_w2 := if (isTRUE(s$use_window2_grid)) ref_minute(sol_all, s$grid_event2) else ref_min_w1]
  sched <- merge(schedules[[p]], refs, by = c("recorder", "date"), all.x = TRUE)
  if (!is.na(row$start1_hour)) sched[, `:=`(w1_start = ref_min_w1 + row$start1_hour * 60,
                                            w1_end = ref_min_w1 + (row$start1_hour + row$dur1_h) * 60)]
  if (!is.na(row$start2_hour)) sched[, `:=`(w2_start = ref_min_w2 + row$start2_hour * 60,
                                            w2_end = ref_min_w2 + (row$start2_hour + row$dur2_h) * 60)]
  keep <- if (row$window_type == "W1") sched[, block_in_window(start_min, end_min, w1_start, w1_end)] else
    if (row$window_type == "W2") sched[, block_in_window(start_min, end_min, w2_start, w2_end)] else
      sched[, block_in_window(start_min, end_min, w1_start, w1_end) | block_in_window(start_min, end_min, w2_start, w2_end)]
  sf <- sched[keep][as.Date(date) %in% unique(as.Date(dt$date))]
  sf[, intersect(c("w1_start", "w1_end", "w2_start", "w2_end", "ref_min_w1", "ref_min_w2"), names(sf)) := NULL]
  res <- compute_richness(dt = copy(dt), schedule = sf, per_spatial = "all", per_temporal = "day", return_species = TRUE)
  all_days <- unique(sf[, .(date = as.Date(date))])
  rm_ <- merge(all_days, as.data.table(res$richness)[, .(date = as.Date(date), richness)], by = "date", all.x = TRUE)
  rm_[is.na(richness), richness := 0]
  c(effort = sum(sf$end_min - sf$start_min + 1L), richness = mean(rm_$richness),
    n_species = length(unique(as.data.table(res$species)$species)))
}

## ---- 5. Plots (same code as SCORE.R, inputs replaced by settings) -----------
fill_scale <- function(pal, name) if (has_scico) scico::scale_fill_scico(palette = pal, name = name, direction = -1) else
  scale_fill_viridis_c(name = name, direction = -1)

plot_raref <- function(tr) {
  ci <- tr[, .(mean_richness = mean(richness_pct, na.rm = TRUE), lower = quantile(richness_pct, 0.025, na.rm = TRUE),
               upper = quantile(richness_pct, 0.975, na.rm = TRUE)), by = effort]
  gg <- ggplot(tr, aes(x = effort, y = richness_pct, group = bootstrap, color = factor(bootstrap))) +
    geom_line(alpha = 0.5, aes(text = paste0("<br>Duty cycle: ", duty_cycle_label, "<br>Effort: ", effort,
                                             " min<br>Richness: ", round(richness_pct, 1), "%<br>Bootstrap ", bootstrap))) +
    geom_point(alpha = 0.5) + labs(x = "Sampling effort [min]", y = "Relative richness [%]") +
    theme_minimal() + guides(color = "none")
  if (isTRUE(s$show_summary))
    gg <- gg + geom_ribbon(data = ci, aes(x = effort, ymin = lower, ymax = upper), inherit.aes = FALSE,
                           fill = "darkblue", alpha = 0.4) +
      geom_line(data = ci, aes(x = effort, y = mean_richness), inherit.aes = FALSE, color = "red")
  gg
}
plot_missed <- function(ms) {
  pmin_ <- min(ms$period); pmax_ <- max(ms$period); pst <- unique(diff(sort(unique(ms$period))))[1]
  if (is.na(pst)) pst <- 1
  ggplot(ms, aes(x = period, y = species_missed, fill = N)) +
    geom_tile(aes(text = paste0("Species: ", species_missed, "\nCycle: ", period, " min\nTimes missed: ", N)),
              color = "white", width = pst, height = 1) +
    scale_fill_viridis_c(option = "D", name = "# bootstraps missed", limits = c(0, max(ms$N, na.rm = TRUE))) +
    scale_x_continuous(name = "Cycle length [min]", breaks = seq(pmin_, pmax_, by = pst), expand = c(0, 0)) +
    labs(y = "Species", title = "Species missed across duty cycles") + theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(angle = 90, vjust = 0.5),
          panel.grid = element_blank())
}
plot_window <- function(df, title, event, pal) {
  df <- copy(df)
  frontier <- df[, .SD[which.max(richness_pct)], by = duration_h]
  df[, tooltip := paste0("<b>", title, "</b><br>Start offset: ", round(start_off_h, 2), " h<br>Duration: ",
                         round(duration_h, 2), " h<br>Relative Richness: ", round(richness_pct, 1),
                         " %<br>Mean daily richness: ", round(richness, 2), "<br>Species: ", n_species_detected)]
  p <- ggplot(df, aes(x = start_off_h, y = duration_h)) +
    geom_tile(aes(fill = richness_pct, text = tooltip), color = "grey89", linewidth = 0.5) +
    geom_contour(aes(z = richness_pct), binwidth = 5, color = "grey70", linewidth = 0.2) +
    fill_scale(pal, "Relative richness [%]") +
    labs(x = paste0("Start offset relative to ", event, " [hours]"), y = "Recording duration [hours]", title = title) +
    theme_minimal() + theme(panel.grid = element_blank(), axis.ticks = element_line()) +
    scale_x_continuous(breaks = seq(floor(min(df$start_off_h)), ceiling(max(df$start_off_h)), by = 0.5)) +
    scale_y_continuous(breaks = seq(floor(min(df$duration_h)), ceiling(max(df$duration_h)), by = 0.5))
  if (isTRUE(s$show_frontier))
    p <- p + geom_point(data = frontier, aes(start_off_h, duration_h), inherit.aes = FALSE,
                        color = "brown1", size = 3, shape = 8)
  p
}
plot_comb <- function(grid, n_recorder_days) {
  df1 <- copy(grid$w1)[, effort_total := duration_h * n_recorder_days * 60]
  df2 <- copy(grid$w2)[, effort_total := duration_h * n_recorder_days * 60]
  df3 <- copy(grid$comb)[, effort_total := effort * n_recorder_days]
  ggplot() +
    geom_point(data = df1, aes(effort_total, richness_pct), color = "steelblue", size = 2.5, alpha = 0.5) +
    geom_point(data = df2, aes(effort_total, richness_pct), color = "orange", size = 2.5, alpha = 0.5) +
    geom_point(data = df3, aes(effort_total, richness_pct), color = "darkgreen", size = 3.5, alpha = 0.7) +
    labs(x = "Total recording effort [min]", y = "Relative richness [%]") + theme_minimal()
}
plot_optimum <- function(dt) {
  pf <- dt[pareto == TRUE][order(effort)]
  dt <- copy(dt)
  dt[, tooltip := paste0("Relative richness: ", round(richness_pct, 1), "%<br>Mean daily richness: ",
                         round(richness, 2), "<br>Type: ", window_type, "<br>Effort: ", effort, " min<br>W1: ",
                         start1_hour, " h / ", dur1_h, " h<br>W2: ", start2_hour, " h / ", dur2_h,
                         " h<br>Period: ", duty_period)]
  ggplot(dt, aes(x = effort, y = richness_pct)) +
    geom_point(aes(color = window_type, shape = window_type, text = tooltip), alpha = 0.7, size = 2) +
    geom_line(data = pf, aes(group = 1), color = "red", linewidth = 0.5) +
    geom_point(data = pf, color = "red", size = 2, alpha = 0.7) +
    scale_color_manual(values = c(W1 = "#1f78b4", W2 = "#6a3d9a", Combined = "#33a02c"), name = "Window type") +
    scale_shape_manual(values = c(W1 = 16, W2 = 17, Combined = 18), name = "Window type") +
    labs(x = "Recording effort [min]", y = "Relative richness [%]",
         title = "Pareto frontier: optimal trade-off richness vs effort") + theme_minimal()
}
save_plot <- function(p, dir, name, w = 8, h = 6) {
  saveRDS(p, file.path(dir, paste0(name, ".rds")))
  suppressWarnings(ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white"))
  if (isTRUE(s$save_plotly) && has_plotly) {
    pl <- suppressWarnings(plotly::ggplotly(p, tooltip = "text"))
    f <- file.path(dir, paste0(name, ".html"))
    ok <- tryCatch({ htmlwidgets::saveWidget(pl, f, selfcontained = TRUE); TRUE }, error = function(e) FALSE)
    if (!ok) htmlwidgets::saveWidget(pl, f, selfcontained = FALSE)
  }
}

## ---- 6. One selection run: Tabs 2-4, verification, saving ----------------------
# Returns the selected schedule; saves everything in folder d and one line in index_file
select_and_save <- function(dt, run_id, d, index_file, extra = list(), return_all = FALSE) {
  t_run <- Sys.time()
  set.seed(s$seed)
  pilot_days <- sort(unique(dt$date))
  rec_dates <- CJ(recorder = unique(dt$recorder), date = c(min(pilot_days) - 1L, pilot_days))
  sol_all <- solar_for(rec_dates)            # pilot days + the day before (schedules crossing midnight)

  t2   <- run_tab2(dt, sol_all);                         msg("  Tab 2 done")
  grid <- run_tab3(dt, sol_all);                         msg("  Tab 3 done")
  cands <- cand_table(window_sets_from(grid))
  res  <- run_tab4(dt, sol_all, t2$schedules, cands);    msg("  Tab 4 done (", nrow(res), " candidates)")

  n_rec <- uniqueN(dt$recorder); n_days <- length(pilot_days)
  ref24 <- mean(dt[, .(n = uniqueN(species)), by = date]$n)
  res[, effort_pct_24h := 100 * effort / (n_days * 1440 * n_rec)]
  res[, effort_reduction_pct := 100 - effort_pct_24h]
  res[, effort_pct_max := 100 * effort / max(effort)]
  res[, richness_pct_24h := 100 * richness / ref24]
  best <- res[richness_pct >= s$threshold][order(effort, -richness_pct)][1]
  near <- res[richness_pct >= s$threshold & effort <= (1 + s$near_optimal_pct / 100) * best$effort][order(effort, -richness_pct)]

  # Verification against the app code path
  vi  <- unique(c(sample.int(nrow(res), min(s$n_verify, nrow(res))),
                  which(res$effort == best$effort & res$richness_pct == best$richness_pct)[1]))
  ver <- rbindlist(lapply(vi, function(i) {
    a <- app_candidate(dt, sol_all, t2$schedules, res[i])
    data.table(row = i, fast_effort = res$effort[i], app_effort = a[["effort"]],
               fast_richness = res$richness[i], app_richness = a[["richness"]],
               fast_n_species = res$n_species_detected[i], app_n_species = a[["n_species"]])
  }))
  max_diff <- ver[, max(abs(fast_effort - app_effort), abs(fast_richness - app_richness),
                        abs(fast_n_species - app_n_species))]
  if (max_diff > 1e-8) { fwrite(ver, file.path(dirname(d), paste0(run_id, "_VERIFICATION_FAILED.csv")))
    stop("Fast engine differs from the app code path for ", run_id, " (max diff ", max_diff, ")") }
  msg("  Verified ", nrow(ver), " candidates against the app code path: identical")

  # ---- Save (same layout as save_run) ----
  dir.create(file.path(d, "plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(d, "tables"), showWarnings = FALSE)
  keep_s <- s[setdiff(names(s), "species_file")]
  params <- data.table(parameter = names(keep_s), value = sapply(keep_s, function(v) paste(v, collapse = ", ")))
  params <- rbind(params, data.table(
    parameter = c("species_file", "data_files", "first_date", "last_date", "n_days", "n_recorders",
                  "recorders", "units", "n_species_in_pilot", "bootstrap_replicates_used",
                  "richness_24h_mean_daily", "n_candidates", "verification_max_diff", "run_date", names(extra)),
    value = c(if (is.null(s$species_file)) "none" else s$species_file, paste(data_files, collapse = "; "),
              format(min(pilot_days)), format(max(pilot_days)), n_days, n_rec,
              paste(unique(dt$recorder), collapse = ", "), paste(sort(unique(dt$unit)), collapse = ", "),
              uniqueN(dt$species), t2$B, round(ref24, 3), nrow(res), max_diff, format(Sys.time()),
              sapply(extra, paste, collapse = ", "))))
  fwrite(params, file.path(d, "parameters.csv"))
  fwrite(best, file.path(d, "best_schedule.csv"))
  fwrite(near, file.path(d, "near_optimal_schedules.csv"))
  fwrite(res,  file.path(d, "tables", "multi_optimum_all_candidates.csv"))
  fwrite(ver,  file.path(d, "tables", "verification_vs_app.csv"))
  saveRDS(grid, file.path(d, "tables", "time_window_data.rds"))
  if (!is.null(t2$summary)) {
    saveRDS(t2$summary, file.path(d, "tables", "duty_cycle_data.rds"))
    fwrite(t2$missed, file.path(d, "tables", "missed_species_summary.csv"))
    save_plot(plot_raref(t2$summary), file.path(d, "plots"), "ggplot_rarefaction")
    save_plot(plot_missed(t2$missed), file.path(d, "plots"), "gg_heatmap_missed", h = max(5, 0.15 * nlevels(t2$missed$species_missed)))
  }
  save_plot(plot_window(grid$w1, "Window 1", s$grid_event, "oslo"), file.path(d, "plots"), "ggplot_w1")
  if (!is.null(grid$w2)) {
    save_plot(plot_window(grid$w2, "Window 2", s$grid_event2, "bilbao"), file.path(d, "plots"), "ggplot_w2")
    save_plot(plot_comb(grid, uniqueN(dt[, .(recorder, date)])), file.path(d, "plots"), "ggplot_comb")
  }
  save_plot(plot_optimum(res), file.path(d, "plots"), "ggplot_optimum")

  idx <- cbind(data.table(run_id = run_id, first_date = min(pilot_days), last_date = max(pilot_days),
                          n_days = n_days, threshold = s$threshold), best)
  idx[, names(idx) := lapply(.SD, as.character)]
  old <- if (file.exists(index_file)) fread(index_file, colClasses = "character") else NULL
  if (!is.null(old)) old <- old[run_id != idx$run_id]          # a re-run replaces its old line
  fwrite(rbind(old, idx, fill = TRUE), index_file)
  msg("  Saved in ", d, " (", round(as.numeric(difftime(Sys.time(), t_run, units = "mins")), 1), " min)")
  if (return_all) return(list(best = best, t2 = t2, grid = grid, res = res, sol_all = sol_all))
  best
}

# Apply ONE schedule to another data set (same app logic: schedules, block_in_window, mean daily richness)
evaluate_schedule <- function(dt, best) {
  days <- sort(unique(dt$date))
  sol  <- solar_for(CJ(recorder = unique(dt$recorder), date = c(min(days) - 1L, days)))
  sch  <- run_tab2(dt, sol, bootstrap = FALSE)$schedules
  sch  <- sch[sapply(sch, attr, "period") == best$duty_period]
  cand <- best[, .(window_type, start1_hour, dur1_h, start2_hour, dur2_h)][, w := 1L]
  r    <- run_tab4(dt, sol, sch, cand)
  a    <- app_candidate(dt, sol, sch, r[1])               # check with the app code path
  if (abs(a[["richness"]] - r$richness[1]) > 1e-8 || abs(a[["effort"]] - r$effort[1]) > 1e-8)
    stop("evaluate_schedule(): fast engine and app code path differ")
  ref24 <- mean(dt[, .(n = uniqueN(species)), by = date]$n)
  data.table(n_days = length(days), n_units = uniqueN(dt$unit), effort = r$effort[1],
             effort_pct_24h = 100 * r$effort[1] / (length(days) * 1440 * uniqueN(dt$recorder)),
             richness = r$richness[1], richness_24h = ref24,
             richness_pct_24h = 100 * r$richness[1] / ref24, n_species = r$n_species_detected[1],
             n_species_24h = uniqueN(dt$species))
}

## ---- 7. CHAPTER A: temporal pilots -------------------------------------------------
if (run_temporal) {
  msg("######## CHAPTER A: temporal pilots ########")
  dir.create(runs_dir, recursive = TRUE, showWarnings = FALSE)
  pilots <- if (!is.null(pilots_custom)) as.data.table(pilots_custom) else
    CJ(start = pilot_starts, length = pilot_lengths, sorted = FALSE)
  pilots[, `:=`(start = as.Date(start), length = as.integer(length))]
  pilots[, end := start + length - 1L]
  pilots[, run_id := sprintf("P%d_%s", length, format(start))]
  msg(nrow(pilots), " pilots: ", paste(pilots$run_id, collapse = ", "))
  for (k in seq_len(nrow(pilots))) {
    pl <- pilots[k]; d <- file.path(runs_dir, pl$run_id)
    if (skip_existing_runs && file.exists(file.path(d, "best_schedule.csv"))) { msg("Skip ", pl$run_id, " (already done)"); next }
    msg("=== ", pl$run_id, " (", format(pl$start), " to ", format(pl$end), ")")
    dt <- dtA[date >= pl$start & date <= pl$end]
    if (nrow(dt) == 0) { msg("No data for ", pl$run_id); next }
    select_and_save(dt, pl$run_id, d, file.path(runs_dir, "index.csv"))
  }
  msg("Chapter A finished. Next: validate_pilots.R with runs_dir = ", runs_dir)
}

## ---- 8. CHAPTER B: spatial validation ---------------------------------------------
if (run_spatial) {
  msg("######## CHAPTER B: spatial validation ########")
  dir.create(spatial_runs_dir, recursive = TRUE, showWarnings = FALSE)
  dS <- dtA[date >= as.Date(spatial_start) & date <= as.Date(spatial_end)]
  units  <- sort(unique(dS$unit))
  groups <- unit_group(units)
  msg(length(units), " recorders with data between ", spatial_start, " and ", spatial_end, ": ",
      paste(units, collapse = ", "))

  # Splits (training units -> the other units are the test set)
  splits <- list()
  if ("group" %in% spatial_splits) for (g in unique(groups))
    splits[[length(splits) + 1]] <- list(type = "group", id = paste0("S_group_", g), train = units[groups == g])
  if ("loo" %in% spatial_splits) for (u in units)
    splits[[length(splits) + 1]] <- list(type = "loo", id = paste0("S_loo_", u), train = setdiff(units, u))
  if ("random" %in% spatial_splits) {
    stopifnot(n_train < length(units))
    set.seed(s$seed); seen <- character(0); i <- 0
    while (i < n_random_splits) {
      tr <- sort(sample(units, n_train)); key <- paste(tr, collapse = "|")
      if (key %in% seen) next
      seen <- c(seen, key); i <- i + 1
      splits[[length(splits) + 1]] <- list(type = "random", id = sprintf("S_random%02d", i), train = tr)
    }
  }
  fwrite(rbindlist(lapply(splits, function(x) data.table(run_id = x$id, split_type = x$type,
         train_units = paste(x$train, collapse = ", "), test_units = paste(setdiff(units, x$train), collapse = ", ")))),
         file.path(spatial_runs_dir, "spatial_splits.csv"))
  msg(length(splits), " splits (listed in spatial_splits.csv)")

  out <- list()
  for (x in splits) {
    test <- setdiff(units, x$train); d <- file.path(spatial_runs_dir, x$id)
    msg("=== ", x$id, ": train ", length(x$train), " recorders, test ", length(test))
    dtr <- dS[unit %in% x$train]; dte <- dS[unit %in% test]
    best <- if (skip_existing_runs && file.exists(file.path(d, "best_schedule.csv"))) {
      msg("  selection already saved, re-using it"); fread(file.path(d, "best_schedule.csv"))
    } else select_and_save(dtr, x$id, d, file.path(spatial_runs_dir, "index.csv"),
                           extra = list(split_type = x$type, train_units = x$train, test_units = test))
    ins <- evaluate_schedule(dtr, best)          # in-sample (training recorders)
    oos <- evaluate_schedule(dte, best)          # out-of-sample (test recorders)
    out[[length(out) + 1]] <- data.table(
      run_id = x$id, split_type = x$type, train_units = paste(x$train, collapse = ", "),
      test_units = paste(test, collapse = ", "), n_train = length(x$train), n_test = length(test),
      window_type = best$window_type, start1_hour = best$start1_hour, dur1_h = best$dur1_h,
      start2_hour = best$start2_hour, dur2_h = best$dur2_h, duty_period = best$duty_period,
      effort_pct_24h = ins$effort_pct_24h,
      in_sample_richness_pct_24h = ins$richness_pct_24h,
      out_of_sample_richness_pct_24h = oos$richness_pct_24h,
      difference = oos$richness_pct_24h - ins$richness_pct_24h,
      in_sample_species = ins$n_species, in_sample_species_24h = ins$n_species_24h,
      out_of_sample_species = oos$n_species, out_of_sample_species_24h = oos$n_species_24h)
    msg("  in-sample ", round(ins$richness_pct_24h, 1), " %  |  out-of-sample ", round(oos$richness_pct_24h, 1), " %")
  }
  sp <- rbindlist(out)
  fwrite(sp, file.path(spatial_runs_dir, "spatial_validation_summary.csv"))
  long <- melt(sp, id.vars = c("run_id", "split_type"),
               measure.vars = c("in_sample_richness_pct_24h", "out_of_sample_richness_pct_24h"),
               variable.name = "set", value.name = "pct")
  long[, set := factor(set, c("in_sample_richness_pct_24h", "out_of_sample_richness_pct_24h"),
                       c("In-sample (training recorders)", "Out-of-sample (test recorders)"))]
  pS <- ggplot(long, aes(set, pct, group = run_id, colour = split_type)) +
    geom_line(alpha = 0.6) + geom_point(size = 2.5) +
    labs(x = NULL, y = "Mean daily richness retained\n(% of continuous 24-h recording)", colour = "Split") +
    theme_minimal()
  ggsave(file.path(spatial_runs_dir, "fig_spatial_in_vs_out.png"), pS, width = 7, height = 5, dpi = 300, bg = "white")
  saveRDS(pS, file.path(spatial_runs_dir, "fig_spatial_in_vs_out.rds"))
  print(sp[, .(run_id, n_train, n_test, effort_pct_24h, in_sample_richness_pct_24h,
               out_of_sample_richness_pct_24h, difference)], digits = 3)
  msg("Chapter B finished. Results: ", file.path(spatial_runs_dir, "spatial_validation_summary.csv"))
}

## ---- 9. CHAPTER C: SCORE presets ------------------------------------------------------
if (run_presets) {
  msg("######## CHAPTER C: SCORE presets ########")
  dir.create(presets_out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a
  s_default <- s
  for (pr in presets_to_make) {
    meta <- modifyList(preset_common, pr[setdiff(names(pr), c("start", "end", "s_override", "validation"))])
    s <- if (is.null(pr$s_override)) s_default else modifyList(s_default, pr$s_override)   # used by run_tab2-4
    dt <- dtA[date >= as.Date(pr$start) & date <= as.Date(pr$end)]
    if (nrow(dt) == 0) { msg("No data between ", pr$start, " and ", pr$end); next }
    dts <- sort(unique(dt$date))
    # preset ID as in the app: country code_region_first date_n days
    code <- if (nzchar(meta$country %||% "")) unname(COUNTRIES[meta$country]) else "XX"
    if (is.na(code)) stop("Country '", meta$country, "' is not in the app list (COUNTRIES)")
    region <- gsub("[^A-Za-z0-9]+", "-", trimws(meta$region %||% ""))
    meta$preset_id <- paste(c(code, if (nzchar(region)) region, format(min(dts)), paste0(length(dts), "d")), collapse = "_")
    zip_file <- file.path(presets_out_dir, paste0("SCORE_preset_", meta$preset_id, ".zip"))
    if (skip_existing_runs && file.exists(zip_file)) { msg("Skip ", meta$preset_id, " (already done)"); next }
    msg("=== Preset ", meta$preset_id, " (", pr$start, " to ", pr$end, ")")

    run <- select_and_save(dt, meta$preset_id, file.path(presets_out_dir, meta$preset_id),
                           file.path(presets_out_dir, "index.csv"), return_all = TRUE)

    settings <- list(
      grid_event = s$grid_event, grid_event2 = s$grid_event2, use_window2_grid = s$use_window2_grid,
      duty_duration = as.numeric(s$duty_duration),
      lat = if (s$coord_mode == "single") s$single_lat else NA, lon = if (s$coord_mode == "single") s$single_lon else NA,
      timezone = s$timezone, min_conf = s$min_conf,
      species_filter = if (is.null(s$species_file)) "none" else paste("species list:", basename(s$species_file)),
      n_units = uniqueN(dt$unit),
      window_grid = list(
        w1 = list(reference_event = s$grid_event, start_offset_h = c(s$grid_start_min, s$grid_start_max),
                  step_h = s$grid_start_step, duration_h = c(s$grid_dur_min, s$grid_dur_max), duration_step_h = s$grid_dur_step),
        w2 = if (isTRUE(s$use_window2_grid)) list(reference_event = s$grid_event2,
                  start_offset_h = c(s$grid_start_min2, s$grid_start_max2), step_h = s$grid_start_step2,
                  duration_h = c(s$grid_dur_min2, s$grid_dur_max2), duration_step_h = s$grid_dur_step2),
        top_n_combined = s$top_n_combined))

    core <- build_preset_core(settings = settings, dt = dt, combo = run$res, solar = run$sol_all,
                              schedules = run$t2$schedules, targets = preset_content$targets,
                              include_pareto = preset_content$include_pareto,
                              include_species = preset_content$include_species,
                              species_min_conf = preset_content$species_min_conf, score_version = SCORE_VERSION)
    preset <- add_preset_meta(core, meta, settings)
    if (!is.null(pr$validation)) preset$validation <- pr$validation
    preset$provenance <- list(generated_by = "dev/SCORE_batch.R (CHAPTER C), same workflow as the app, verified against the app code path",
                              data_files = basename(data_files))

    # check: the 90 % schedule of the preset is the best schedule of the run (same rule)
    r90 <- Filter(function(z) z$target_pct == s$threshold, preset$results$recommended_schedules)
    if (length(r90) && isTRUE(r90[[1]]$reached) && abs(r90[[1]]$mean_daily_richness - round(run$best$richness, 2)) > 0.01)
      stop("Preset and run disagree on the ", s$threshold, " % schedule for ", meta$preset_id)

    n_rd <- uniqueN(dt[, .(recorder, date)])
    plots <- list(
      duty_cycle_rarefaction    = if (!is.null(run$t2$summary)) plot_raref(run$t2$summary),
      duty_cycle_missed_species = if (!is.null(run$t2$missed))  plot_missed(run$t2$missed),
      time_window_w1            = plot_window(run$grid$w1, "Window 1", s$grid_event, "oslo"),
      time_window_w2            = if (!is.null(run$grid$w2)) plot_window(run$grid$w2, "Window 2", s$grid_event2, "bilbao"),
      time_window_combined      = if (!is.null(run$grid$w2)) plot_comb(run$grid, n_rd),
      multi_optimum             = plot_optimum(run$res))
    write_score_preset_zip(preset, plots, zip_file)
    if (!is.null(library_presets_dir)) {
      utils::unzip(zip_file, exdir = library_presets_dir)           # dossier <id>/ avec .json + figures/
      msg("  Preset unzipped into ", library_presets_dir)
    }
    msg("  Preset saved: ", zip_file)
  }
  s <- s_default
  msg("Chapter C finished. Unzip the presets into preset_library/presets/ and run render_all_presets().")
}
msg("Finished in ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min.")
