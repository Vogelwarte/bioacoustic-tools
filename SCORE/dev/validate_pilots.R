# =============================================================================
# validate_pilots.R
# Out-of-sample (temporal) validation of SCORE pilot schedules
#
# For every pilot run saved with save_run() (one folder per run containing
# best_schedule.csv + parameters.csv), the script:
#   1. rebuilds the selected schedule exactly as the Multi-optimum tab does
#      (duty-cycle blocks over 24 h, kept when they fall inside W1 / W2,
#       using block_in_window() from functions.R);
#   2. checks that it reproduces the app value on the pilot days;
#   3. evaluates it on the other days of the season, relative to continuous
#      24-h recording on the same days:
#        - mean daily richness retained (%), total species retained (%),
#        - per block of N days (seasonal transferability),
#        - species missed (and whether they were already present in the pilot);
#   4. saves all numbers, settings and figures in <runs_dir>/validation/.
#
# Window anchoring (column "anchor" in the outputs):
#   "solar" = windows follow sunrise/sunset of each day (what the ms describes;
#             MAIN RESULT)
#   "fixed" = windows fixed at the clock times of the pilot (mean pilot
#             sunrise/sunset), i.e. a recorder that cannot follow solar times
#   "app"   = windows anchored to the sunrise/sunset of the FIRST DAY OF THE
#             DATA LOADED IN THE APP. This is what the original app does: the
#             original block_in_window() only uses the first row of the
#             schedule, and the schedule covers every date of the imported
#             data (not only the pilot days). Used for the reproduction check
#             of runs made with the original functions.R.
#
# Mode: detections of all recorders are pooled ("global" mode, as in the
# runs). Only the prepared detection table is cached; everything else is
# recomputed each time (it takes seconds).
# =============================================================================

## ---- 0. USER SETTINGS (edit here) -------------------------------------------
cfg <- list(
  runs_dir          = "/Users/amandine.serrurier/Documents/SCORE_runs/batch_runs",
  out_dir           = "/Users/amandine.serrurier/Documents/SCORE_runs/validation",
  functions_R       = "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/VoWa/shiny_App_raref/bioacoustic-tools/SCORE/functions.R",
  data_file         = "/Users/amandine.serrurier/Documents/traning_data_complete_bn_2024_GG/PVP_2023_24h/GON_GRE_2023.csv",
  species_file      = NULL,            # path to "Scientific_Common" list, or NULL = no filter (as in your runs)
  exclude_recorders = c("RUB6-ROLLE"),
  min_conf          = 0,                 # same as the app runs
  lat               = 46.19,           # SAME single coordinates as in the app
  lon               = 8.13,
  tz                = "Europe/Zurich",   # same time zone as in the app
  season_start      = "2023-05-12",
  season_end        = "2023-09-14",
  block_days        = 14,                # length of the seasonal test blocks
  date_from         = "file",            # "file" = date of the file start (as the app);
                                         # "detection" = real date of the detection
  # run_pattern       = "^1205-1805_7day$",  # regex to select run folders (NULL = all)
  anchors           = c("solar", "fixed", "app"),
  app_version       = "fixed",             # "old"   = runs made with the ORIGINAL block_in_window()
                                         # "fixed" = runs made with the corrected block_in_window()
  app_first_date    = NULL               # first date of the data imported in the app
                                         # (NULL = first date in data_file)
)
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(suncalc)
  library(ggplot2)
})
stopifnot(!is.na(cfg$lat), !is.na(cfg$lon), file.exists(cfg$functions_R), file.exists(cfg$data_file))
source(cfg$functions_R)                 # create_sun_duty_schedule(), block_in_window()
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time()
msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)

## ---- 1. Detections: prepare once, then cache --------------------------------
cache_file <- file.path(cfg$out_dir, "prepared_detections.rds")
cache_key  <- list(
  data = file.info(cfg$data_file)[, c("size", "mtime")],
  species = if (!is.null(cfg$species_file)) file.info(cfg$species_file)[, c("size", "mtime")] else NULL,
  cfg[c("exclude_recorders", "min_conf", "date_from", "species_file")],
  version = "prep_v2"            # change when the preparation code changes
)

prepare_detections <- function() {
  msg("Reading ", basename(cfg$data_file), " ...")
  hdr  <- names(fread(cfg$data_file, nrows = 0))
  # Same column and same parsing as the app (SCORE.R, dt_parsed): timestamp taken
  # from the FIRST "YYYYMMDD_HHMMSS" in basename(Begin Path) without extension
  fcol <- if ("Begin Path" %in% hdr) "Begin Path" else "Begin File"
  x <- fread(cfg$data_file, select = c(fcol, "Begin Time (s)", "Common Name", "Confidence"),
             showProgress = FALSE)
  setnames(x, c("file", "offset_s", "species", "conf"))
  n_raw <- nrow(x)

  # Diagnostic: do Begin File and Begin Path give the same timestamp?
  if (all(c("Begin File", "Begin Path") %in% hdr)) {
    bf <- fread(cfg$data_file, select = "Begin File", showProgress = FALSE)[[1]]
    first_ts <- function(v) regmatches(v, regexpr("\\d{8}_\\d{6}", v))[match(seq_along(v), which(grepl("\\d{8}_\\d{6}", v)))]
    ub <- unique(data.table(bf = basename(bf), bp = x$file))
    ub[, `:=`(t_bf = first_ts(bf), t_bp = first_ts(tools::file_path_sans_ext(basename(gsub("\\\\+", "/", bp)))))]
    n_diff <- ub[is.na(t_bf) != is.na(t_bp) | (t_bf != t_bp), .N]
    msg("Files where Begin File and Begin Path give different timestamps: ", n_diff, " of ", nrow(ub))
    rm(bf, ub)
  }

  x[, file := tools::file_path_sans_ext(basename(gsub("\\\\+", "/", file)))]   # as normalize_begin_path() + basename
  x[, recorder := sub("_\\d{8}_\\d{6}.*$", "", file)]
  uf  <- unique(x$file)
  hit <- regexpr("\\d{8}_\\d{6}", uf)                                   # first match, as str_extract()
  uts <- ifelse(hit > 0, substr(uf, hit, hit + 14L), NA_character_)
  x[, ts := uts[match(file, uf)]]
  x <- x[!is.na(ts) & !(recorder %in% cfg$exclude_recorders) & conf >= cfg$min_conf]

  if (!is.null(cfg$species_file)) {
    sp <- trimws(sub("\r$", "", readLines(cfg$species_file, encoding = "UTF-8", warn = FALSE)))
    sp <- unique(sub("^[^_]*_", "", sp[sp != ""]))
    not_found <- setdiff(sp, unique(x$species))
    if (length(not_found)) msg("Species in list never detected: ", paste(not_found, collapse = ", "))
    x <- x[species %in% sp]
  }

  # Integer arithmetic instead of POSIXct (fast; no DST change between May and Sep)
  u <- unique(x[, .(ts)])
  u[, `:=`(fdate = as.Date(substr(ts, 1, 8), "%Y%m%d"),
           sod   = as.integer(substr(ts, 10, 11)) * 3600L +
                   as.integer(substr(ts, 12, 13)) * 60L + as.integer(substr(ts, 14, 15)))]
  x <- u[x, on = "ts"]
  x[, tot_s := sod + floor(as.numeric(offset_s))]
  x[, minute := as.integer((tot_s %/% 60) %% 1440)]          # = hour*60 + minute of detection
  x[, date := if (cfg$date_from == "file") fdate else fdate + as.integer(tot_s %/% 86400)]
  n_shift <- x[fdate != fdate + as.integer(tot_s %/% 86400), .N]

  species_levels <- sort(unique(x$species))
  D <- unique(x[, .(date, minute, sp = match(species, species_levels))])
  setkey(D, date, minute)
  list(D = D, species = species_levels, recorders = sort(unique(x$recorder)),
       n_raw = n_raw, n_kept = nrow(x), n_after_midnight = n_shift)
}

if (file.exists(cache_file) && identical(readRDS(cache_file)$key, cache_key)) {
  msg("Using cached detections")
  prep <- readRDS(cache_file)$prep
} else {
  prep <- prepare_detections()
  saveRDS(list(key = cache_key, prep = prep), cache_file)
}
D <- prep$D[date >= as.Date(cfg$season_start) & date <= as.Date(cfg$season_end)]
season_days <- sort(unique(D$date))
msg(format(nrow(D), big.mark = ","), " unique date-minute-species records, ",
    length(season_days), " days, ", length(prep$species), " species")
if (prep$n_after_midnight > 0)
  msg("NOTE: ", prep$n_after_midnight, " detections lie after midnight of their file start date ",
      "(date_from = '", cfg$date_from, "')")

# Continuous 24-h reference: species per day
REF_day <- D[, .(ref = uniqueN(sp)), keyby = date]

## ---- 2. Solar table (same function and rounding as the app) -----------------
app_first <- if (is.null(cfg$app_first_date)) min(prep$D$date) else as.Date(cfg$app_first_date)
sun_days  <- seq(min(app_first, min(season_days)), max(season_days), by = "day")
sun <- as.data.table(getSunlightTimes(date = sun_days, lat = cfg$lat, lon = cfg$lon,
                                      keep = c("dawn", "sunrise", "sunset", "dusk"), tz = cfg$tz))
sun[, date := as.Date(date)]
for (ev in c("dawn", "sunrise", "sunset", "dusk")) {
  sun[[paste0(ev, "_min")]] <- as.integer(hour(sun[[ev]]) * 60 + minute(sun[[ev]]))
}

ref_minutes <- function(event, days) {
  if (event %in% c("dawn", "sunrise", "sunset", "dusk")) {
    sun[data.table(date = days), on = "date"][[paste0(event, "_min")]]
  } else {                                           # fixed clock event "HH:MM"
    rep(as.integer(sub(":.*", "", event)) * 60L + as.integer(sub(".*:", "", event)), length(days))
  }
}

## ---- 3. Read the pilot runs -------------------------------------------------
run_dirs <- list.dirs(cfg$runs_dir, recursive = FALSE)
run_dirs <- run_dirs[file.exists(file.path(run_dirs, "best_schedule.csv")) &
                     file.exists(file.path(run_dirs, "parameters.csv"))]
if (!is.null(cfg$run_pattern)) run_dirs <- run_dirs[grepl(cfg$run_pattern, basename(run_dirs))]
stopifnot(length(run_dirs) > 0)

runs <- rbindlist(lapply(run_dirs, function(d) {
  p <- fread(file.path(d, "parameters.csv"), colClasses = "character")
  p <- setNames(as.list(p$value), p$parameter)
  b <- fread(file.path(d, "best_schedule.csv"))
  num <- function(v) if (is.null(v) || length(v) == 0 || is.na(v) || v == "") NA_real_ else as.numeric(v)
  data.table(run_id = basename(d),
             first = as.Date(p$first_date), last = as.Date(p$last_date),
             n_days = as.integer(p$n_days),
             window_type = b$window_type,
             ev1 = p$grid_event, off1 = num(b$start1_hour), dur1 = num(b$dur1_h),
             ev2 = p$grid_event2, off2 = num(b$start2_hour), dur2 = num(b$dur2_h),
             period = as.numeric(b$duty_period), duty = as.numeric(p$duty_duration),
             app_effort = as.numeric(b$effort), app_richness = as.numeric(b$richness),
             app_richness_pct = as.numeric(b$richness_pct),
             params_signature = paste(p[setdiff(names(p), c("first_date", "last_date", "n_days"))],
                                      collapse = "|"))
}))
setorder(runs, first, n_days)
if (uniqueN(runs$params_signature) > 1)
  warning("Runs do not all use the same settings (see parameters.csv of each run).")
msg(nrow(runs), " pilot runs found")

## ---- 4. Build and evaluate a schedule ---------------------------------------
# Recorded minutes for a run on a set of days, with a given window anchoring
build_schedule <- function(r, days, anchor) {
  # Duty-cycle blocks over 24 h, exactly as create_sun_duty_schedule(record_24h = TRUE)
  base <- create_sun_duty_schedule(sun_table = sun[1, .(date, dawn, sunrise, sunset, dusk)],
                                   period = r$period, duty_duration = r$duty,
                                   record_24h = TRUE)
  blocks <- unique(base[, .(start_min, end_min)])
  S <- CJ(date = days, b = seq_len(nrow(blocks)))
  S[, `:=`(start_min = blocks$start_min[b], end_min = blocks$end_min[b])]

  pilot_days <- seq(r$first, r$last, by = "day")
  ref_for <- function(ev) {
    switch(anchor,
           solar = ref_minutes(ev, S$date),
           fixed = rep(as.integer(round(mean(ref_minutes(ev, pilot_days)))), nrow(S)),
           app   = rep(ref_minutes(ev, app_first), nrow(S)))
  }
  use1 <- r$window_type %in% c("W1", "Combined")
  use2 <- r$window_type %in% c("W2", "Combined")
  if (use1) { rr <- ref_for(r$ev1); S[, `:=`(w1s = rr + r$off1 * 60, w1e = rr + (r$off1 + r$dur1) * 60)] }
  if (use2) { rr <- ref_for(r$ev2); S[, `:=`(w2s = rr + r$off2 * 60, w2e = rr + (r$off2 + r$dur2) * 60)] }

  # block_in_window() uses scalar window bounds -> apply per day
  S[, keep := {
    k <- rep(FALSE, .N)
    if (use1) k <- k | block_in_window(start_min, end_min, w1s[1], w1e[1])
    if (use2) k <- k | block_in_window(start_min, end_min, w2s[1], w2e[1])
    k
  }, by = date]
  S <- S[keep == TRUE]
  effort <- S[, .(effort = sum(end_min - start_min + 1)), keyby = date]
  mins <- if (all(S$end_min == S$start_min)) S[, .(date, minute = start_min)] else
    S[, .(minute = start_min:end_min), by = .(date, start_min)][, .(date, minute)]
  list(minutes = unique(mins), effort = effort)
}

# Metrics on a set of days
evaluate <- function(sch, days) {
  if (length(days) == 0) return(NULL)
  Dd  <- D[date %in% days]
  hit <- Dd[sch$minutes[date %in% days], on = .(date, minute), nomatch = 0L]
  dd  <- data.table(date = days)
  rd  <- hit[, .(rich = uniqueN(sp)), keyby = date][dd, on = "date"]
  rd[is.na(rich), rich := 0L]
  ref <- REF_day[dd, on = "date"][is.na(ref), ref := 0L]
  eff <- sch$effort[dd, on = "date"][is.na(effort), effort := 0]
  sp_sched <- unique(hit$sp); sp_ref <- unique(Dd$sp)
  list(
    metrics = data.table(
      n_days            = length(days),
      mean_rich_sched   = mean(rd$rich),
      mean_rich_24h     = mean(ref$ref),
      daily_richness_pct_24h = 100 * mean(rd$rich) / mean(ref$ref),
      n_species_sched   = length(sp_sched),
      n_species_24h     = length(sp_ref),
      total_species_pct_24h = 100 * length(sp_sched) / length(sp_ref),
      effort_pct_24h    = 100 * sum(eff$effort) / (length(days) * 1440)),
    sp_sched = sp_sched, sp_ref = sp_ref)
}

## ---- 5. Loop over runs and anchors ------------------------------------------
summary_rows <- list(); block_rows <- list(); missed_rows <- list(); check_rows <- list()
check_anchor <- if (cfg$app_version == "old") "app" else "solar"
if (!check_anchor %in% cfg$anchors) cfg$anchors <- c(cfg$anchors, check_anchor)
msg("Reproduction check uses anchor '", check_anchor, "' (app_version = '", cfg$app_version,
    "'; app data assumed to start on ", format(app_first), ")")
blocks_of <- function(days) as.integer(as.numeric(days - as.Date(cfg$season_start)) %/% cfg$block_days) + 1L

for (i in seq_len(nrow(runs))) {
  r <- runs[i]
  pilot_days <- season_days[season_days >= r$first & season_days <= r$last]
  after_days <- season_days[season_days > r$last]
  other_days <- season_days[season_days < r$first | season_days > r$last]

  for (anc in cfg$anchors) {
    sch <- build_schedule(r, season_days, anc)
    ins <- evaluate(sch, pilot_days)

    if (anc == check_anchor) {               # reproduction check vs the app
      check_rows[[length(check_rows) + 1]] <- data.table(
        run_id = r$run_id, app_richness = r$app_richness,
        script_richness = ins$metrics$mean_rich_sched,
        app_effort = r$app_effort,
        script_effort = sum(sch$effort[data.table(date = pilot_days), on = "date"]$effort, na.rm = TRUE))
    }

    for (set_name in c("pilot_in_sample", "after_pilot", "all_other_days")) {
      days <- switch(set_name, pilot_in_sample = pilot_days, after_pilot = after_days,
                     all_other_days = other_days)
      ev <- if (set_name == "pilot_in_sample") ins else evaluate(sch, days)
      if (is.null(ev)) next
      summary_rows[[length(summary_rows) + 1]] <- cbind(
        r[, .(run_id, pilot_start = first, pilot_end = last, pilot_length = n_days, window_type,
              off1, dur1, off2, dur2, period)],
        anchor = anc, day_set = set_name, ev$metrics)
    }

    # Seasonal blocks (pilot days excluded)
    for (bk in unique(blocks_of(other_days))) {
      days <- other_days[blocks_of(other_days) == bk]
      ev <- evaluate(sch, days)
      block_rows[[length(block_rows) + 1]] <- cbind(
        r[, .(run_id, pilot_start = first, pilot_length = n_days)], anchor = anc,
        block = bk, block_start = min(days), block_end = max(days),
        position = ifelse(min(days) > r$last, "after_pilot", "before_pilot"),
        weeks_after_pilot_end = ifelse(min(days) > r$last,
                                       as.numeric(mean(days) - r$last) / 7, NA_real_),
        ev$metrics)
    }

    # Species missed after the pilot (main anchor only)
    if (anc == "solar" && length(after_days)) {
      ev <- evaluate(sch, after_days)
      missed <- setdiff(ev$sp_ref, ev$sp_sched)
      in_pilot <- unique(D[date %in% pilot_days]$sp)
      if (length(missed))
        missed_rows[[length(missed_rows) + 1]] <- data.table(
          run_id = r$run_id, pilot_start = r$first, pilot_length = r$n_days,
          species = prep$species[missed],
          detected_in_pilot_24h = missed %in% in_pilot,
          n_days_detected_24h_after = sapply(missed, function(s) uniqueN(D[date %in% after_days & sp == s]$date)))
    }
  }
  msg("Run ", i, "/", nrow(runs), " done: ", r$run_id)
}

summary_dt <- rbindlist(summary_rows)
blocks_dt  <- rbindlist(block_rows)
missed_dt  <- rbindlist(missed_rows)
check_dt   <- rbindlist(check_rows)
check_dt[, `:=`(richness_diff = script_richness - app_richness, effort_diff = script_effort - app_effort)]

## ---- 6. Save numbers and settings -------------------------------------------
fwrite(summary_dt, file.path(cfg$out_dir, "validation_summary.csv"))
fwrite(blocks_dt,  file.path(cfg$out_dir, "validation_by_block.csv"))
fwrite(missed_dt,  file.path(cfg$out_dir, "validation_missed_species.csv"))
fwrite(check_dt,   file.path(cfg$out_dir, "check_reproduction_vs_app.csv"))
settings <- data.table(setting = names(cfg), value = sapply(cfg, function(v) paste(v, collapse = ", ")))
settings <- rbind(settings, data.table(
  setting = c("run_date", "n_runs", "runs", "n_detections_raw", "n_detections_kept",
              "n_detections_after_midnight", "n_species", "recorders", "n_season_days", "R_version"),
  value = c(format(Sys.time()), nrow(runs), paste(runs$run_id, collapse = ", "), prep$n_raw,
            prep$n_kept, prep$n_after_midnight, length(prep$species),
            paste(prep$recorders, collapse = ", "), length(season_days), R.version.string)))
fwrite(settings, file.path(cfg$out_dir, "validation_settings.csv"))

## ---- 7. Figures -------------------------------------------------------------
save_fig <- function(p, name, w = 8, h = 5) {
  ggsave(file.path(cfg$out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  saveRDS(p, file.path(cfg$out_dir, paste0(name, ".rds")))
}
main <- summary_dt[anchor == "solar" & day_set %in% c("pilot_in_sample", "after_pilot")]
main[, day_set := factor(day_set, c("pilot_in_sample", "after_pilot"),
                         c("In-sample (pilot days)", "Out-of-sample (later days)"))]
lab_levels <- format(sort(unique(runs$first)), "%d %b")
start_lab  <- function(d) factor(format(d, "%d %b"), levels = lab_levels)   # facets in date order
main[, pilot_start := start_lab(pilot_start)]

pA <- ggplot(main, aes(pilot_length, daily_richness_pct_24h, colour = day_set)) +
  geom_line() + geom_point(size = 2.5) +
  scale_x_continuous(breaks = sort(unique(main$pilot_length))) +
  facet_wrap(~ pilot_start, nrow = 1) +
  labs(x = "Pilot length (days)", y = "Mean daily richness retained\n(% of continuous 24-h recording)",
       colour = NULL) + theme_minimal() + theme(legend.position = "bottom")
save_fig(pA, "fig_A_in_vs_out_of_sample", w = 10)

pB <- ggplot(blocks_dt[anchor == "solar" & position == "after_pilot"],
             aes(weeks_after_pilot_end, daily_richness_pct_24h, colour = factor(pilot_length))) +
  geom_line() + geom_point() +
  facet_wrap(~ start_lab(pilot_start), nrow = 1) +
  labs(x = "Weeks after the end of the pilot", y = "Mean daily richness retained (% of 24 h)",
       colour = "Pilot length (d)") + theme_minimal()
save_fig(pB, "fig_B_transferability_over_season", w = 10)

pC <- ggplot(summary_dt[day_set == "after_pilot"],
             aes(factor(pilot_length), daily_richness_pct_24h, colour = anchor)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  facet_wrap(~ start_lab(pilot_start), nrow = 1) +
  labs(x = "Pilot length (days)", y = "Out-of-sample richness retained (% of 24 h)",
       colour = "Windows") + theme_minimal()
save_fig(pC, "fig_C_solar_vs_fixed_windows", w = 10)

if (nrow(missed_dt)) {
  mt <- missed_dt[, .(n_runs = .N, share_absent_from_pilot = mean(!detected_in_pilot_24h)),
                  by = .(species, pilot_length)]
  pD <- ggplot(mt, aes(factor(pilot_length), reorder(species, n_runs), fill = n_runs)) +
    geom_tile() + labs(x = "Pilot length (days)", y = NULL, fill = "Runs where\nspecies missed") +
    theme_minimal()
  save_fig(pD, "fig_D_missed_species", w = 7, h = max(4, 0.18 * uniqueN(mt$species)))
}

## ---- 8. Console summary -----------------------------------------------------
cat("\n=== Reproduction check (should be ~0) ===\n"); print(check_dt)
cat("\n=== Main results (solar windows) ===\n")
main_tab <- dcast(summary_dt[anchor == "solar"], run_id + pilot_start + pilot_length ~ day_set,
                  value.var = "daily_richness_pct_24h")
main_tab <- merge(main_tab, summary_dt[anchor == "solar" & day_set == "pilot_in_sample",
                                       .(run_id, effort_pct_24h)], by = "run_id")
setcolorder(main_tab, c("run_id", "pilot_start", "pilot_length", "effort_pct_24h", "pilot_in_sample"))
print(main_tab[order(pilot_start, pilot_length)], digits = 3)
msg("Finished in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s. Outputs in ", cfg$out_dir)

