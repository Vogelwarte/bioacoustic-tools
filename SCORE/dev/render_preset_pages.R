# =============================================================================
# render_preset_pages.R
# Turns SCORE presets into readable GitHub pages (Markdown + figures) and a
# catalogue (presets/README.md).
#
# Usage:
#   setwd("<folder that contains the presets/ folder>")
#   source("render_preset_pages.R")
#   render_all_presets("presets")
#
# Input: presets/ can contain
#   - unzipped presets downloaded from the app:  presets/<id>/<id>.json + presets/<id>/figures/*.png
#   - or single .json files:                      presets/<anything>.json
# Output:
#   presets/README.md                       catalogue (table of all presets)
#   presets/index.csv                       same catalogue, machine-readable
#   presets/pages/<id>.md                   one readable page per preset
#   presets/pages/<id>_schedule.png         recommended schedules, relative to sunrise/sunset
#   presets/pages/<id>_curve.png            effort-richness curve
#
# Packages: jsonlite, ggplot2 (jsonlite is called as jsonlite::fromJSON() and is
# not attached on purpose: library(jsonlite) masks shiny::validate()).
# =============================================================================
setwd("/Users/amandine.serrurier/Documents/SCORE_runs/github_presets")
source("render_preset_pages.R")
render_all_presets("presets")
library(ggplot2)

# ---- small helpers -----------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Missing, empty or "TO_FILL" values are shown as an en dash
val <- function(x, suffix = "") {
  if (is.null(x) || length(x) == 0) return("–")
  x <- unlist(x)
  x <- trimws(as.character(x[!is.na(x)]))
  x <- x[x != "" & !grepl("^TO_FILL", x)]
  if (!length(x)) return("–")
  paste0(paste(x, collapse = ", "), suffix)
}

fmt_h <- function(h) {                      # 1.5 -> "1 h 30", 2 -> "2 h"
  hh <- floor(h); mm <- round((h - hh) * 60)
  paste0(hh, " h", if (mm > 0) sprintf(" %02d", mm) else "")
}

fmt_offset <- function(event, offset_h) {
  if (offset_h == 0) return(event)
  paste0(event, if (offset_h > 0) " + " else " − ", fmt_h(abs(offset_h)))
}

fmt_window <- function(w) {
  start <- fmt_offset(w$reference_event, w$start_offset_h)
  end   <- fmt_offset(w$reference_event, w$start_offset_h + w$duration_h)
  paste0("from ", start, " to ", end, " (≈ ", sub("-", "–", val(w$approx_clock_time)), ")")
}

fmt_duty <- function(label) {
  if (identical(label, "continuous")) return("continuous")
  x <- as.numeric(strsplit(sub(" min", "", label), "/")[[1]])
  paste0(x[1], " min every ", x[2], " min")
}

# hour bin k (k to k+1 h after sunrise) as text
fmt_bin <- function(k) {
  if (k >= 0) paste0(k, "–", k + 1, " h after sunrise") else
    paste0(abs(k + 1), "–", abs(k), " h before sunrise")
}

# text sparkline of the activity profile (24 bins: 6 h before to 18 h after sunrise)
sparkline <- function(share) {
  if (!length(share) || sum(share) == 0) return("–")
  ch <- c("▁", "▂", "▃", "▄", "▅", "▆", "▇", "█")
  s <- ifelse(share == 0, "·", ch[pmax(1, ceiling(8 * share / max(share)))])
  paste0("`", paste(s[1:6], collapse = ""), "┃", paste(s[7:24], collapse = ""), "`")
}

md_table <- function(df) {
  df[] <- lapply(df, function(x) gsub("\\|", "\\\\|", as.character(x)))
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}

details <- function(summary, body, open = FALSE)
  c(if (open) "<details open>" else "<details>", paste0("<summary><b>", summary, "</b></summary>"), "",
    body, "", "</details>", "")

n_recorders <- function(p) p$sampling$recorders$n %||% p$sampling$recorders$n_units %||% NA

reached <- function(p) Filter(function(z) !isFALSE(z$reached) && length(z$windows),
                              p$results$recommended_schedules)

# ---- figures -----------------------------------------------------------------

# Recommended schedules, in hours relative to the reference event (sunrise / sunset):
# the bars do not move with the season, unlike clock times.
plot_schedule <- function(p, file) {
  recs <- reached(p)
  rows <- list()
  for (z in recs) {
    lab <- paste0(z$target_pct, " %  (", fmt_duty(z$duty_cycle$label), ")")
    for (w in z$windows) rows[[length(rows) + 1]] <- data.frame(
      target = lab, event = w$reference_event,
      start = w$start_offset_h, end = w$start_offset_h + w$duration_h)
  }
  d <- do.call(rbind, rows)
  d$target <- factor(d$target, levels = rev(unique(d$target)))
  d$panel  <- factor(paste("Hours relative to", d$event), levels = unique(paste("Hours relative to", d$event)))
  zero <- unique(d[, c("panel", "event")])

  g <- ggplot(d) +
    geom_vline(data = zero, aes(xintercept = 0), linetype = "dashed", colour = "#e0a400") +
    geom_text(data = zero, aes(x = 0, y = Inf, label = event), vjust = -0.5, size = 3.3, colour = "#9a7000") +
    geom_segment(aes(x = start, xend = end, y = target, yend = target), linewidth = 6, colour = "#1f6f8b") +
    facet_wrap(~panel, scales = "free_x") +
    scale_x_continuous(breaks = function(l) seq(floor(l[1]), ceiling(l[2]), by = 1),
                       labels = function(x) ifelse(x > 0, paste0("+", x), as.character(x))) +
    coord_cartesian(clip = "off") +
    labs(x = "Hours", y = NULL, title = "When to record",
         subtitle = "Recording windows of the recommended schedules (same data as the table)") +
    theme_minimal(base_size = 11) +
    theme(plot.background = element_rect(fill = "white", colour = NA),
          panel.grid.minor = element_blank(), panel.spacing = unit(1.5, "lines"),
          strip.text = element_text(face = "bold", margin = margin(b = 14)))
  ggsave(file, g, width = 8, height = 1.9 + 0.42 * nlevels(d$target), dpi = 150, bg = "white")
}

# Effort-richness curve (Pareto front of all candidate schedules tested)
plot_curve <- function(p, file) {
  pf <- do.call(rbind, lapply(p$results$pareto_front, function(z)
    data.frame(min_day = z$effort_pct_24h * 14.4, richness = z$richness_pct_of_24h)))
  pf <- pf[pf$min_day > 0, ]
  recs <- reached(p)
  rec <- do.call(rbind, lapply(recs, function(z)
    data.frame(min_day = z$effort_pct_24h * 14.4, richness = z$richness_pct_of_24h,
               label = paste0(z$target_pct, " %"))))
  best <- if (length(recs)) recs[[1]]$richness_pct_of_24h / recs[[1]]$richness_pct_of_best * 100 else NA
  brk <- c(1, 5, 15, 30, 60, 120, 240, 480, 1440)

  g <- ggplot(pf, aes(min_day, richness)) +
    geom_hline(yintercept = 100, colour = "grey30", linewidth = 0.4) +
    annotate("text", x = min(pf$min_day), y = 100, label = "continuous 24 h recording",
             hjust = 0, vjust = -0.5, size = 3.2, colour = "grey30") +
    geom_step(colour = "grey55", linewidth = 0.6) +
    geom_point(colour = "grey55", size = 0.9) +
    geom_point(data = rec, colour = "#1f6f8b", size = 2.8) +
    geom_text(data = rec, aes(label = label), colour = "#1f6f8b", hjust = 1.3, vjust = -0.5, size = 3.4) +
    scale_x_log10(breaks = brk, labels = function(x) ifelse(x >= 60, paste0(x / 60, " h"), paste0(x, " min"))) +
    scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25), labels = function(y) paste0(y, " %")) +
    labs(x = "Recording time per day and recorder (log scale)",
         y = "Mean daily species richness\n(% of 24 h recording)",
         subtitle = "Grey: best schedule found for each recording time. Blue: recommended schedules.") +
    theme_minimal(base_size = 11) +
    theme(plot.background = element_rect(fill = "white", colour = NA), panel.grid.minor = element_blank())
  if (!is.na(best) && best < 99) g <- g +
    geom_hline(yintercept = best, linetype = "dashed", colour = "#1f6f8b", linewidth = 0.4) +
    annotate("text", x = min(pf$min_day), y = best, label = "best schedule tested",
             hjust = 0, vjust = -0.5, size = 3.2, colour = "#1f6f8b")
  ggsave(file, g, width = 7.5, height = 4.4, dpi = 150, bg = "white")
}

# ---- one page ----------------------------------------------------------------

render_preset_page <- function(json_file, presets_dir) {
  p  <- jsonlite::fromJSON(json_file, simplifyVector = FALSE)
  id <- p$preset_id
  pages_dir <- file.path(presets_dir, "pages")
  dir.create(pages_dir, showWarnings = FALSE, recursive = TRUE)
  message("Rendering ", id)
  rel <- function(f) paste0("../", sub(paste0("^", normalizePath(presets_dir), "/?"), "", normalizePath(f)))

  plot_schedule(p, file.path(pages_dir, paste0(id, "_schedule.png")))
  plot_curve(p,    file.path(pages_dir, paste0(id, "_curve.png")))

  st <- p$study; sa <- p$sampling; cl <- p$classifier; an <- p$analysis
  n_days <- sa$period$n_days
  recs <- reached(p)
  targets <- vapply(recs, function(z) z$target_pct, 0)

  # --- header ---------------------------------------------------------------
  authors <- val(vapply(p$contributors, function(a) a$name %||% "", ""))
  head <- c(
    paste0("# ", trimws(val(p$title))), "",
    paste0("`", id, "` · ", authors, " · SCORE ", val(p$score_version)), "",
    paste0("[⬇ Download the preset (.json, for the SCORE app)](", rel(json_file), ") · ",
           "[← Back to catalogue](../README.md)"), "",
    paste0("> **Pilot data:** ", n_days, " days (", sa$period$first_date, " to ", sa$period$last_date, "), ",
           n_recorders(p), " recorder(s), ", val(an$reference$total_species %||% an$reference$total_species_24h),
           " species detected with continuous 24 h recording (", round(an$reference$continuous_24h_mean_daily_richness, 1),
           " species per day on average)."), "")

  # --- recommended schedules --------------------------------------------------
  rec <- do.call(rbind, lapply(p$results$recommended_schedules, function(z) {
    if (isFALSE(z$reached)) return(data.frame(`Target` = paste0(z$target_pct, " %"), `Duty cycle` = "not reached",
      `Time windows` = "–", `Recording per day` = "–", `Richness (% of 24 h)` = "–",
      `Richness (% of best tested)` = "–", `Species never detected` = "–", check.names = FALSE))
    om <- z$omission
    data.frame(
      `Target` = paste0(z$target_pct, " %"),
      `Duty cycle` = fmt_duty(z$duty_cycle$label),
      `Time windows` = paste(vapply(z$windows, fmt_window, ""), collapse = "<br>"),
      `Recording per day` = paste0(round(z$recording_h_per_day * 60), " min (", z$effort_pct_24h, " % of 24 h)"),
      `Richness (% of 24 h)` = paste0(z$richness_pct_of_24h, " %"),
      `Richness (% of best tested)` = paste0(z$richness_pct_of_best, " %"),
      `Species never detected` = if (is.null(om)) "–" else
        paste0(om$species_never_detected, " of ", om$species_detected_24h),
      check.names = FALSE)
  }))
  best24 <- if (length(recs)) round(recs[[1]]$richness_pct_of_24h / recs[[1]]$richness_pct_of_best * 100, 1) else NA

  sec_rec <- c(
    "## Recommended schedules", "",
    paste0("For each target, the schedule that needs the **least recording** among the ",
           format(an$n_candidates, big.mark = " "), " schedules tested by SCORE."), "",
    md_table(rec), "",
    details("How to read this table", c(
      "- **Target**: share of species richness to reach, relative to the *best schedule tested* (see below).",
      "- **Duty cycle**: how to record *inside* the time windows. *Continuous* = record without pause; *1 min every 7 min* = record 1 min, pause 6 min. It is the duty cycle of the least-effort schedule found, not a fixed choice: when a short continuous window is cheaper than a long sampled one, SCORE recommends the continuous window.",
      "- **Time windows**: relative to sunrise or sunset, so that they follow the season. Clock times (≈) are only indicative: they use the mean sunrise/sunset of the pilot days.",
      "- **Recording per day**: per recorder.",
      "- **Richness (% of 24 h)**: mean number of species detected per day with this schedule, as a % of the species detected per day with continuous 24 h recording (all recorders pooled per day).",
      paste0("- **Richness (% of best tested)**: same, but relative to the best schedule among those tested. SCORE only tests windows within the grid chosen in the app (Tab 3), so the best schedule tested can detect fewer species than 24 h recording", if (!is.na(best24)) paste0(" (here: ", best24, " % of 24 h)") else "", ". The targets refer to this value."),
      "- **Species never detected**: species detected at least once during the pilot with 24 h recording but never with this schedule (details below).")),
    paste0("![When to record](", id, "_schedule.png)"), "",
    paste0("*Each bar is a recording window of a recommended schedule (same data as the table). Positions are in hours relative to ",
           "sunrise or sunset, so the windows follow the change of day length; the duty cycle in brackets applies inside the bars.*"), "")

  # --- omission -------------------------------------------------------------
  oms <- Filter(function(z) !is.null(z$omission), recs)
  sec_om <- c("## Species missed by each schedule", "")
  if (length(oms)) {
    om_tab <- do.call(rbind, lapply(oms, function(z) data.frame(
      `Target` = paste0(z$target_pct, " %"),
      `Recording per day` = paste0(round(z$recording_h_per_day * 60), " min"),
      `Species detected over the pilot` = paste0(z$omission$species_detected, " of ", z$omission$species_detected_24h),
      `Species never detected` = paste0(z$omission$species_never_detected, " (", z$omission$species_never_detected_pct, " %)"),
      `Mean daily omission` = paste0(z$omission$mean_daily_omission_pct, " %"),
      check.names = FALSE)))
    lists <- unlist(lapply(oms, function(z) details(
      paste0(z$target_pct, " % schedule: ", z$omission$species_never_detected, " species never detected"),
      if (length(z$omission$missed_species)) paste(unlist(z$omission$missed_species), collapse = ", ") else "_none_")))
    sec_om <- c(sec_om, md_table(om_tab), "",
      "- **Species never detected**: detected at least once with 24 h recording during the pilot, but on no day with this schedule.",
      "- **Mean daily omission**: on an average day, share of the species detected with 24 h recording that the schedule misses.",
      "- Detections per species and per schedule are in the table *Species* below.", "", lists)
  } else sec_om <- c(sec_om, "_Not available in this preset (created with an earlier version of SCORE)._", "")

  # --- species --------------------------------------------------------------
  sec_sp <- character(0)
  if (length(p$species)) {
    spd <- do.call(rbind, lapply(p$species, function(s) {
      share <- unlist(s$activity_share_by_hour_from_sunrise %||% s$activity_profile_h_from_sunrise$share)
      peak  <- if (length(share) && sum(share) > 0) fmt_bin((-6:17)[which.max(share)]) else "–"
      by_s  <- s$days_detected_by_schedule %||%
        if (!is.null(s$days_detected_recommended_90)) list(`90%` = s$days_detected_recommended_90) else list()
      cols <- setNames(lapply(paste0(targets, "%"), function(k) val(by_s[[k]])), paste0(targets, " %"))
      data.frame(
        Species = s$common_name,
        Detections = val(s$n_detections %||% s$n_detection_minutes),
        `Mean confidence` = val(s$confidence$mean),
        `24 h` = s$days_detected_24h,
        as.data.frame(cols, check.names = FALSE),
        `Peak activity` = peak,
        `Activity over the day` = sparkline(share),
        check.names = FALSE)
    }))
    if (all(spd$`Mean confidence` == "–")) spd$`Mean confidence` <- NULL
    spd <- spd[order(-spd$`24 h`, spd$Species), ]
    det_note <- if (is.null(p$species[[1]]$n_detections)) "minutes with at least one detection" else "number of BirdNET detections"
    sec_sp <- c(
      "## Species", "",
      paste0("All species detected during the pilot (", nrow(spd), "), with ", val(cl$name),
             if (val(cl$min_confidence) != "–") paste0(" (confidence threshold: ", cl$min_confidence, ")") else "",
             ". Species names are those of the classifier; detections were not verified unless stated in the notes."), "",
      details(paste0("Show the table of the ", nrow(spd), " species"), c(
        paste0("- **Detections**: ", det_note, " during the pilot, all recorders."),
        "- **Mean confidence**: mean classifier confidence of these detections (low values point to possible false positives).",
        paste0("- **24 h, ", paste(paste0(targets, " %"), collapse = ", "), "**: number of pilot days (out of ", n_days,
               ") on which the species was detected with continuous 24 h recording and with each recommended schedule."),
        "- **Peak activity**: the hour, relative to sunrise, with the most detections.",
        "- **Activity over the day**: detections per hour from 6 h before sunrise to 18 h after it; `┃` marks sunrise, `·` no detection.", "",
        md_table(spd)), open = TRUE))
  }

  # --- curve ------------------------------------------------------------------
  sec_curve <- c(
    "## How many species for how much recording?", "",
    paste0("![Effort-richness curve](", id, "_curve.png)"), "",
    paste0("*Each grey point is the best schedule found by SCORE for a given recording time per day (among ",
           format(an$n_candidates, big.mark = " "), " schedules tested). The curve shows how fast richness increases with recording time ",
           "and where adding recording time brings few additional species. Use it to choose a compromise other than the recommended targets.*"), "")

  # --- context --------------------------------------------------------------
  sites <- if (length(st$sites)) paste(vapply(st$sites, function(s) paste0(s$name, " (", s$n_recorders, ")"), ""), collapse = ", ") else NULL
  habitat <- val(c(unlist(st$habitat$codes %||% st$habitat$code), st$habitat$description))
  elev <- if (!is.null(st$elevation_mean_m)) val(st$elevation_mean_m, " m (mean)") else
    if (val(st$elevation_m$min) == "–") "–" else paste0(st$elevation_m$min, "–", st$elevation_m$max, " m")
  loc <- st$location %||% st$location_rounded
  loc_txt <- if (is.null(loc)) "not shared" else
    paste0(loc$lat, "° N, ", loc$lon, "° E", if (isTRUE(loc$rounded) || !is.null(loc$precision_deg)) " (rounded to ~10 km)" else "")
  ctx <- data.frame(
    Item = c("Project", "Objective", "Target taxa", "Target species", "Country / region", "Sites (recorders)",
             "Habitat (EUNIS)", "Elevation", "Location", "Period", "Timing", "Recorders",
             "Sample rate", "Deployment height", "Classifier", "Confidence threshold"),
    Value = c(
      val(st$project_description %||% st$project), val(st$objective), val(st$target_taxa), val(st$target_species),
      paste(val(st$country), "/", val(st$region)), val(sites), habitat, elev, loc_txt,
      paste0(sa$period$first_date, " to ", sa$period$last_date, " (", n_days, " days with data, ", val(sa$period$timezone), ")"),
      val(sa$period$timing %||% sa$period$season),
      paste0(n_recorders(p), " × ", val(sa$recorders$model)),
      val(sa$recorders$sample_rate_hz, " Hz"), val(sa$recorders$deployment_height_m, " m"),
      paste(val(cl$name), val(cl$version)), val(cl$min_confidence)))
  ctx <- ctx[ctx$Value != "–", ]
  sec_ctx <- c("## Study context", "", md_table(ctx), "",
               if (!is.null(cl$note)) c(paste0("_Note: ", cl$note, "_"), "") else NULL,
               if (val(p$notes) != "–") c(paste0("_Notes: ", p$notes, "_"), "") else NULL)

  # --- validation (optional) --------------------------------------------------
  v <- p$validation
  val_lines <- character(0)
  if (!is.null(v$temporal)) val_lines <- c(val_lines,
    "**Temporal transfer** (90 % schedule applied to other days, same recorders):", "",
    md_table(data.frame(
      Test = c("Pilot days (in-sample)", "Later days of the season", "All other days of the season"),
      `Richness (% of 24 h)` = paste0(c(val(v$temporal$in_sample_pct_of_24h), val(v$temporal$later_days_pct_of_24h),
                                         val(v$temporal$all_other_days_pct_of_24h)), " %"),
      check.names = FALSE)), "")
  if (!is.null(v$spatial)) {
    sp <- v$spatial
    rng <- function(x) if (is.null(x)) "–" else sprintf("%.1f %% (%.1f–%.1f)", x$mean, x$min, x$max)
    val_lines <- c(val_lines,
      "**Spatial transfer** (schedule selected on some recorders, applied to others, same days):", "",
      md_table(data.frame(
        Test = c("Between sites (mean of both directions)", "Random recorder subsets", "Single recorder left out"),
        `Richness (% of 24 h)` = c(
          if (is.null(sp$between_sites)) "–" else paste0(round(mean(unlist(sp$between_sites)), 1), " %"),
          rng(sp$random_10_to_4), rng(sp$leave_one_recorder_out)),
        check.names = FALSE)), "")
  }
  if (!is.null(v$guidance)) val_lines <- c(val_lines, paste0("> ", v$guidance), "")
  if (!length(val_lines)) val_lines <- "_No validation provided with this preset._"
  sec_val <- c("## How well does it transfer?", "", val_lines, "")

  # --- method -----------------------------------------------------------------
  dc <- an$duty_cycles
  cyc <- unlist(dc$cycle_min %||% c(dc$cycle_min_from, dc$cycle_min_to))
  sec_met <- c("## Method", "",
    paste0("- Richness: ", val(an$richness_metric), "; ", val(an$spatial_pooling)),
    paste0("- Duty cycles tested: ", dc$record_min, " min every ", min(cyc), " to ", max(cyc), " min"),
    paste0("- Schedules tested: ", format(an$n_candidates, big.mark = " "),
           " (duty cycles × time windows of the grid of Tab 3)"),
    paste0("- Selection rule: ", val(an$selection_rule)), "")

  # --- figures of the app -----------------------------------------------------
  sec_fig <- character(0)
  figs <- list.files(file.path(dirname(json_file), "figures"), pattern = "\\.png$", full.names = TRUE)
  if (length(figs)) {
    titles <- c(
      duty_cycle_rarefaction = "Tab 2 – Richness by duty cycle (bootstrap)",
      duty_cycle_missed_species = "Tab 2 – Species missed by duty cycle",
      time_window_w1 = "Tab 3 – Richness by time window (window 1)",
      time_window_w2 = "Tab 3 – Richness by time window (window 2)",
      time_window_combined = "Tab 3 – Richness of combined windows",
      multi_optimum = "Tab 4 – Multi-optimum: effort vs richness of all schedules tested")
    key <- sub("\\.png$", "", basename(figs))
    figs <- figs[order(match(key, names(titles)))]; key <- sub("\\.png$", "", basename(figs))
    body <- unlist(lapply(seq_along(figs), function(i) c(
      paste0("**", titles[key[i]] %||% key[i], "**"), "", paste0("![", key[i], "](", rel(figs[i]), ")"), "")))
    sec_fig <- c("## Figures from the SCORE app", "",
                 details("Show the original figures of this analysis", body))
  }

  sec_cite <- c("## Cite", "", val(p$citation %||%
    "Serrurier et al. (in revision) SCORE: a user-friendly application to optimise recording schedules in passive acoustic monitoring. Journal of Avian Biology."), "")

  md <- c(head, sec_rec, sec_om, sec_sp, sec_curve, sec_ctx, sec_val, sec_met, sec_fig, sec_cite)
  writeLines(md, file.path(pages_dir, paste0(id, ".md")), useBytes = TRUE)

  r90 <- Filter(function(z) z$target_pct == 90, recs)
  data.frame(
    preset_id = id, title = trimws(val(p$title)), country = val(st$country),
    habitat = habitat, timing = val(sa$period$timing %||% sa$period$season),
    first_date = sa$period$first_date, n_days = n_days, n_recorders = n_recorders(p),
    min_day_90 = if (length(r90)) round(r90[[1]]$recording_h_per_day * 60) else NA,
    duty_90 = if (length(r90)) fmt_duty(r90[[1]]$duty_cycle$label) else NA,
    stringsAsFactors = FALSE)
}

# ---- catalogue ---------------------------------------------------------------

render_all_presets <- function(presets_dir = "presets") {
  files <- list.files(presets_dir, pattern = "\\.json$", full.names = TRUE, recursive = TRUE)
  files <- files[!grepl("(^|/)pages/", files)]
  if (!length(files)) stop("No .json preset in ", presets_dir)
  idx <- do.call(rbind, lapply(files, render_preset_page, presets_dir = presets_dir))
  idx <- idx[order(idx$country, idx$first_date), ]
  write.csv(idx, file.path(presets_dir, "index.csv"), row.names = FALSE)

  tab <- data.frame(
    Preset = paste0("[", idx$title, "](pages/", idx$preset_id, ".md)"),
    Country = idx$country, Habitat = idx$habitat, Timing = idx$timing,
    Start = idx$first_date, Days = idx$n_days, Recorders = idx$n_recorders,
    `90 % target reached with` = ifelse(is.na(idx$min_day_90), "–",
                                       paste0(idx$min_day_90, " min/day (", idx$duty_90, ")")),
    check.names = FALSE)

  md <- c(
    "# SCORE preset library", "",
    "Recording schedules optimised with [SCORE](../README.md) and shared by users.",
    "Click a preset to see its recommended schedules, the species missed by each schedule, the study context and validation.", "",
    paste0("**", nrow(idx), " preset(s)** · last update ", Sys.Date()), "",
    md_table(tab), "",
    "## Share your own preset", "",
    "1. In SCORE, open the tab **6. Share results**, fill in the form and download the preset (.zip).",
    "2. Open a new issue with the template **Submit a SCORE preset** and attach the .zip, or send it by e-mail.",
    "3. After a quick check, it is added here and archived with a DOI at each release.", "")
  writeLines(md, file.path(presets_dir, "README.md"), useBytes = TRUE)
  message("Done: ", nrow(idx), " preset(s) -> ", normalizePath(presets_dir))
  invisible(idx)
}
