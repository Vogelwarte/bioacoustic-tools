# =============================================================================
# preset_functions.R  --  SCORE presets (shareable summary of an analysis)
# Paste at the END of functions.R (or source() it after functions.R).
# Needs: data.table, lubridate, jsonlite
# =============================================================================

SCORE_PRESET_SCHEMA <- "0.2"

# ---- Vocabularies for the form ----------------------------------------------

OBJECTIVES <- c("community (species richness)" = "community",
                "target species" = "target species")      # users can also type their own

TAXA <- c("birds", "bats", "amphibians", "orthopterans", "mammals")   # users can also type their own

TIMING <- list(
  Season = c("breeding", "post-breeding", "autumn migration", "wintering",
             "spring migration", "year-round"),
  Month  = month.name)

# Country names (ISO 3166-1, from the ISOcodes package) -> ISO alpha-2 code
COUNTRIES <- c(
  "Afghanistan" = "AF", "\u00c5land Islands" = "AX", "Albania" = "AL", "Algeria" = "DZ",
  "American Samoa" = "AS", "Andorra" = "AD", "Angola" = "AO", "Anguilla" = "AI",
  "Antarctica" = "AQ", "Antigua and Barbuda" = "AG", "Argentina" = "AR", "Armenia" = "AM",
  "Aruba" = "AW", "Australia" = "AU", "Austria" = "AT", "Azerbaijan" = "AZ",
  "Bahamas" = "BS", "Bahrain" = "BH", "Bangladesh" = "BD", "Barbados" = "BB",
  "Belarus" = "BY", "Belgium" = "BE", "Belize" = "BZ", "Benin" = "BJ",
  "Bermuda" = "BM", "Bhutan" = "BT", "Bolivia" = "BO", "Bonaire, Sint Eustatius and Saba" = "BQ",
  "Bosnia and Herzegovina" = "BA", "Botswana" = "BW", "Bouvet Island" = "BV", "Brazil" = "BR",
  "British Indian Ocean Territory" = "IO", "Brunei Darussalam" = "BN", "Bulgaria" = "BG", "Burkina Faso" = "BF",
  "Burundi" = "BI", "Cabo Verde" = "CV", "Cambodia" = "KH", "Cameroon" = "CM",
  "Canada" = "CA", "Cayman Islands" = "KY", "Central African Republic" = "CF", "Chad" = "TD",
  "Chile" = "CL", "China" = "CN", "Christmas Island" = "CX", "Cocos (Keeling) Islands" = "CC",
  "Colombia" = "CO", "Comoros" = "KM", "Congo" = "CG", "Congo, The Democratic Republic of the" = "CD",
  "Cook Islands" = "CK", "Costa Rica" = "CR", "C\u00f4te d\'Ivoire" = "CI", "Croatia" = "HR",
  "Cuba" = "CU", "Cura\u00e7ao" = "CW", "Cyprus" = "CY", "Czechia" = "CZ",
  "Denmark" = "DK", "Djibouti" = "DJ", "Dominica" = "DM", "Dominican Republic" = "DO",
  "Ecuador" = "EC", "Egypt" = "EG", "El Salvador" = "SV", "Equatorial Guinea" = "GQ",
  "Eritrea" = "ER", "Estonia" = "EE", "Eswatini" = "SZ", "Ethiopia" = "ET",
  "Falkland Islands (Malvinas)" = "FK", "Faroe Islands" = "FO", "Fiji" = "FJ", "Finland" = "FI",
  "France" = "FR", "French Guiana" = "GF", "French Polynesia" = "PF", "French Southern Territories" = "TF",
  "Gabon" = "GA", "Gambia" = "GM", "Georgia" = "GE", "Germany" = "DE",
  "Ghana" = "GH", "Gibraltar" = "GI", "Greece" = "GR", "Greenland" = "GL",
  "Grenada" = "GD", "Guadeloupe" = "GP", "Guam" = "GU", "Guatemala" = "GT",
  "Guernsey" = "GG", "Guinea" = "GN", "Guinea-Bissau" = "GW", "Guyana" = "GY",
  "Haiti" = "HT", "Heard Island and McDonald Islands" = "HM", "Holy See (Vatican City State)" = "VA", "Honduras" = "HN",
  "Hong Kong" = "HK", "Hungary" = "HU", "Iceland" = "IS", "India" = "IN",
  "Indonesia" = "ID", "Iran" = "IR", "Iraq" = "IQ", "Ireland" = "IE",
  "Isle of Man" = "IM", "Israel" = "IL", "Italy" = "IT", "Jamaica" = "JM",
  "Japan" = "JP", "Jersey" = "JE", "Jordan" = "JO", "Kazakhstan" = "KZ",
  "Kenya" = "KE", "Kiribati" = "KI", "Kuwait" = "KW", "Kyrgyzstan" = "KG",
  "Laos" = "LA", "Latvia" = "LV", "Lebanon" = "LB", "Lesotho" = "LS",
  "Liberia" = "LR", "Libya" = "LY", "Liechtenstein" = "LI", "Lithuania" = "LT",
  "Luxembourg" = "LU", "Macao" = "MO", "Madagascar" = "MG", "Malawi" = "MW",
  "Malaysia" = "MY", "Maldives" = "MV", "Mali" = "ML", "Malta" = "MT",
  "Marshall Islands" = "MH", "Martinique" = "MQ", "Mauritania" = "MR", "Mauritius" = "MU",
  "Mayotte" = "YT", "Mexico" = "MX", "Micronesia, Federated States of" = "FM", "Moldova" = "MD",
  "Monaco" = "MC", "Mongolia" = "MN", "Montenegro" = "ME", "Montserrat" = "MS",
  "Morocco" = "MA", "Mozambique" = "MZ", "Myanmar" = "MM", "Namibia" = "NA",
  "Nauru" = "NR", "Nepal" = "NP", "Netherlands" = "NL", "New Caledonia" = "NC",
  "New Zealand" = "NZ", "Nicaragua" = "NI", "Niger" = "NE", "Nigeria" = "NG",
  "Niue" = "NU", "Norfolk Island" = "NF", "North Korea" = "KP", "North Macedonia" = "MK",
  "Northern Mariana Islands" = "MP", "Norway" = "NO", "Oman" = "OM", "Pakistan" = "PK",
  "Palau" = "PW", "Palestine, State of" = "PS", "Panama" = "PA", "Papua New Guinea" = "PG",
  "Paraguay" = "PY", "Peru" = "PE", "Philippines" = "PH", "Pitcairn" = "PN",
  "Poland" = "PL", "Portugal" = "PT", "Puerto Rico" = "PR", "Qatar" = "QA",
  "R\u00e9union" = "RE", "Romania" = "RO", "Russian Federation" = "RU", "Rwanda" = "RW",
  "Saint Barth\u00e9lemy" = "BL", "Saint Helena, Ascension and Tristan da Cunha" = "SH", "Saint Kitts and Nevis" = "KN", "Saint Lucia" = "LC",
  "Saint Martin (French part)" = "MF", "Saint Pierre and Miquelon" = "PM", "Saint Vincent and the Grenadines" = "VC", "Samoa" = "WS",
  "San Marino" = "SM", "Sao Tome and Principe" = "ST", "Saudi Arabia" = "SA", "Senegal" = "SN",
  "Serbia" = "RS", "Seychelles" = "SC", "Sierra Leone" = "SL", "Singapore" = "SG",
  "Sint Maarten (Dutch part)" = "SX", "Slovakia" = "SK", "Slovenia" = "SI", "Solomon Islands" = "SB",
  "Somalia" = "SO", "South Africa" = "ZA", "South Georgia and the South Sandwich Islands" = "GS", "South Korea" = "KR",
  "South Sudan" = "SS", "Spain" = "ES", "Sri Lanka" = "LK", "Sudan" = "SD",
  "Suriname" = "SR", "Svalbard and Jan Mayen" = "SJ", "Sweden" = "SE", "Switzerland" = "CH",
  "Syria" = "SY", "Taiwan" = "TW", "Tajikistan" = "TJ", "Tanzania" = "TZ",
  "Thailand" = "TH", "Timor-Leste" = "TL", "Togo" = "TG", "Tokelau" = "TK",
  "Tonga" = "TO", "Trinidad and Tobago" = "TT", "Tunisia" = "TN", "T\u00fcrkiye" = "TR",
  "Turkmenistan" = "TM", "Turks and Caicos Islands" = "TC", "Tuvalu" = "TV", "Uganda" = "UG",
  "Ukraine" = "UA", "United Arab Emirates" = "AE", "United Kingdom" = "GB", "United States" = "US",
  "United States Minor Outlying Islands" = "UM", "Uruguay" = "UY", "Uzbekistan" = "UZ", "Vanuatu" = "VU",
  "Venezuela" = "VE", "Vietnam" = "VN", "Virgin Islands, British" = "VG", "Virgin Islands, U.S." = "VI",
  "Wallis and Futuna" = "WF", "Western Sahara" = "EH", "Yemen" = "YE", "Zambia" = "ZM",
  "Zimbabwe" = "ZW"
)


# Habitat vocabulary (EUNIS habitat classification, level 1)
EUNIS_L1 <- c(
  "A - Marine habitats"                                   = "A",
  "B - Coastal habitats"                                  = "B",
  "C - Inland surface waters"                             = "C",
  "D - Mires, bogs and fens"                              = "D",
  "E - Grasslands and lands dominated by forbs, mosses or lichens" = "E",
  "F - Heathland, scrub and tundra"                       = "F",
  "G - Woodland, forest and other wooded land"            = "G",
  "H - Inland unvegetated or sparsely vegetated habitats" = "H",
  "I - Regularly or recently cultivated agricultural, horticultural and domestic habitats" = "I",
  "J - Constructed, industrial and other artificial habitats" = "J",
  "X - Habitat complexes"                                 = "X"
)

minutes_to_hm <- function(m) sprintf("%02d:%02d", (round(m) %/% 60) %% 24, round(m) %% 60)

round_coord <- function(x, precision) {
  if (is.null(precision) || is.na(precision)) return(NULL)
  round(round(x / precision) * precision, 6)
}

# Least-effort schedule reaching a richness target (% of best candidate); ties: highest richness
best_schedule_at <- function(combo, target) {
  x <- as.data.table(combo)[richness_pct >= target][order(effort, -richness_pct)]
  if (nrow(x) == 0) NULL else x[1]
}

# Blocks of one candidate schedule, exactly as the Multi-optimum tab keeps them
# schedule: data.table(recorder, date, start_min, end_min); refs: data.table(recorder, date, ref_min_w1, ref_min_w2)
candidate_blocks <- function(schedule, refs, row, dates) {
  s <- merge(as.data.table(schedule), refs, by = c("recorder", "date"), all.x = TRUE)
  s <- s[as.Date(date) %in% as.Date(dates)]
  k <- rep(FALSE, nrow(s))
  if (!is.na(row$start1_hour)) k <- k | block_in_window(s$start_min, s$end_min,
        s$ref_min_w1 + row$start1_hour * 60, s$ref_min_w1 + (row$start1_hour + row$dur1_h) * 60)
  if (!is.na(row$start2_hour)) k <- k | block_in_window(s$start_min, s$end_min,
        s$ref_min_w2 + row$start2_hour * 60, s$ref_min_w2 + (row$start2_hour + row$dur2_h) * 60)
  s[k, .(recorder, date = as.Date(date), start_min, end_min)]
}

# One schedule row -> list for the JSON
schedule_to_list <- function(r, n_days, n_rec, ref24, ev1, ev2, mean_ref1, mean_ref2, duty_duration) {
  eff_pct <- 100 * r$effort / (n_days * 1440 * n_rec)
  w <- list()
  if (!is.na(r$start1_hour)) w[[length(w) + 1]] <- list(
    reference_event = ev1, start_offset_h = r$start1_hour, duration_h = r$dur1_h,
    approx_clock_time = paste0(minutes_to_hm(mean_ref1 + 60 * r$start1_hour), "-",
                               minutes_to_hm(mean_ref1 + 60 * (r$start1_hour + r$dur1_h))))
  if (!is.na(r$start2_hour)) w[[length(w) + 1]] <- list(
    reference_event = ev2, start_offset_h = r$start2_hour, duration_h = r$dur2_h,
    approx_clock_time = paste0(minutes_to_hm(mean_ref2 + 60 * r$start2_hour), "-",
                               minutes_to_hm(mean_ref2 + 60 * (r$start2_hour + r$dur2_h))))
  list(duty_cycle = list(record_min = duty_duration, cycle_min = r$duty_period,
                         label = if (r$duty_period == duty_duration) "continuous" else paste0(duty_duration, "/", r$duty_period, " min")),
       windows = w,
       recording_h_per_day = round(eff_pct * 24 / 100, 2),
       effort_pct_24h      = round(eff_pct, 2),
       richness_pct_of_best = round(r$richness_pct, 1),
       richness_pct_of_24h  = round(100 * r$richness / ref24, 1),
       mean_daily_richness  = round(r$richness, 2))
}

# -----------------------------------------------------------------------------
# Build the preset (a nested list, written to JSON by write_score_preset)
#   meta      : list of form fields (see SCORE.R, Tab "Share results")
#   settings  : list of analysis settings (inputs of Tabs 1-4)
#   dt        : working detections (dt_raw): species, conf, recorder, date, minute_of_day
#   combo     : Multi-optimum table (multi_combo())
#   solar     : solar table (recorder, date, dawn, sunrise, sunset, dusk)
#   schedules : Tab 2 schedules (values$schedules)
# -----------------------------------------------------------------------------
# The preset is built in two steps, so that the app does not recompute the results each time
# a form field changes: build_preset_core() (results, slow) + add_preset_meta() (form, instant).
build_score_preset <- function(meta, settings, dt, combo, solar, schedules,
                               targets = c(75, 80, 85, 90, 95),
                               include_pareto = TRUE, include_species = TRUE,
                               species_min_conf = 0, score_version = "unknown") {
  core <- build_preset_core(settings, dt, combo, solar, schedules, targets, include_pareto,
                            include_species, species_min_conf, score_version)
  add_preset_meta(core, meta, settings)
}

build_preset_core <- function(settings, dt, combo, solar, schedules,
                              targets = c(75, 80, 85, 90, 95),
                              include_pareto = TRUE, include_species = TRUE,
                              species_min_conf = 0, score_version = "unknown") {
  meta <- list()   # form fields are added later by add_preset_meta()
  dt <- as.data.table(dt); dt[, date := as.Date(date)]
  combo <- as.data.table(combo)
  dates <- sort(unique(dt$date)); n_days <- length(dates); n_rec <- uniqueN(dt$recorder)
  ref24 <- mean(dt[, .(n = uniqueN(species)), by = date]$n)

  sol <- as.data.table(solar)[as.Date(date) %in% dates]
  ev_min <- function(ev) if (ev %in% c("dawn", "sunrise", "sunset", "dusk"))
    as.integer(hour(sol[[ev]]) * 60 + minute(sol[[ev]])) else
    rep(as.integer(sub(":.*", "", ev)) * 60L + as.integer(sub(".*:", "", ev)), nrow(sol))
  ev1 <- settings$grid_event
  ev2 <- if (isTRUE(settings$use_window2_grid)) settings$grid_event2 else settings$grid_event
  refs <- sol[, .(recorder, date = as.Date(date))]
  refs[, `:=`(ref_min_w1 = ev_min(ev1), ref_min_w2 = ev_min(ev2), sunrise_min = ev_min("sunrise"))]
  mean_ref1 <- mean(refs$ref_min_w1); mean_ref2 <- mean(refs$ref_min_w2)
  dd <- as.numeric(settings$duty_duration)
  to_list <- function(r) schedule_to_list(r, n_days, n_rec, ref24, ev1, ev2, mean_ref1, mean_ref2, dd)

  periods <- sapply(schedules, attr, "period")
  refs_w <- refs[, .(recorder, date, ref_min_w1, ref_min_w2)]

  # (recorder, date, species) detected by one candidate schedule (row of combo)
  # det: unique detection minutes data.table(recorder, date, m, species)
  detected_by <- function(r, det) {
    b <- candidate_blocks(schedules[[match(r$duty_period, periods)]], refs_w, r, dates)
    b[, `:=`(bs = as.numeric(start_min), be = as.numeric(end_min))]
    unique(b[det, on = .(recorder, date, bs <= m, be >= m), nomatch = 0L,
             .(date = i.date, species = i.species)])
  }
  det_all <- unique(dt[, .(recorder, date, m = as.numeric(minute_of_day), species)])
  sp_day_24 <- unique(det_all[, .(date, species)])

  # Recommended schedules + species omission (all detections, as for richness)
  rec_rows <- lapply(targets, function(t) best_schedule_at(combo, t))
  recs <- lapply(seq_along(targets), function(i) {
    r <- rec_rows[[i]]
    if (is.null(r)) return(list(target_pct = targets[i], reached = FALSE))
    hit <- detected_by(r, det_all)
    missed <- sort(setdiff(unique(sp_day_24$species), unique(hit$species)))
    daily <- merge(sp_day_24[, .(n24 = .N), by = date], hit[, .(n = .N), by = date], by = "date", all.x = TRUE)
    daily[is.na(n), n := 0]
    c(list(target_pct = targets[i], reached = TRUE), to_list(r),
      list(omission = list(
        species_detected_24h = uniqueN(sp_day_24$species),
        species_detected     = uniqueN(hit$species),
        species_never_detected = length(missed),
        species_never_detected_pct = round(100 * length(missed) / uniqueN(sp_day_24$species), 1),
        mean_daily_omission_pct = round(100 * mean(1 - daily$n / daily$n24), 1),
        missed_species = as.list(missed))))
  })

  # Pareto front (dominance as in the app plot)
  pareto <- NULL
  if (include_pareto) {
    o <- combo[order(effort)]
    best_prev <- c(-Inf, cummax(o$richness_pct))[match(o$effort, o$effort)]
    same_max  <- ave(o$richness_pct, o$effort, FUN = max)
    pf <- o[!(best_prev >= richness_pct | same_max > richness_pct)]
    pareto <- lapply(seq_len(nrow(pf)), function(i) { r <- pf[i]; list(
      effort_pct_24h = round(100 * r$effort / (n_days * 1440 * n_rec), 2),
      richness_pct_of_best = round(r$richness_pct, 1), richness_pct_of_24h = round(100 * r$richness / ref24, 1),
      window_type = r$window_type, cycle_min = r$duty_period,
      w1 = if (!is.na(r$start1_hour)) c(r$start1_hour, r$dur1_h) else NULL,
      w2 = if (!is.na(r$start2_hour)) c(r$start2_hour, r$dur2_h) else NULL) })
  }

  # Species-level summaries (detections with confidence >= species_min_conf)
  species <- NULL
  if (include_species) {
    d_raw <- dt[conf >= species_min_conf]
    det_min <- unique(d_raw[, .(recorder, date, m = as.numeric(minute_of_day), species)])
    # detection under each duty cycle alone (24-h blocks of the Tab 2 schedules)
    p_day <- rbindlist(lapply(seq_along(schedules), function(i) {
      s <- as.data.table(schedules[[i]])[as.Date(date) %in% dates, .(recorder, date = as.Date(date),
                                                                   bs = as.numeric(start_min), be = as.numeric(end_min))]
      h <- s[det_min, on = .(recorder, date, bs <= m, be >= m), nomatch = 0L, .(date = i.date, species = i.species)]
      unique(h)[, .(p = uniqueN(date) / n_days), by = species][, period := periods[i]]
    }))
    # days detected with each recommended schedule
    by_sched <- rbindlist(lapply(seq_along(targets), function(i) {
      r <- rec_rows[[i]]; if (is.null(r)) return(NULL)
      detected_by(r, det_min)[, .(days = uniqueN(date)), by = species][, target := targets[i]]
    }))
    # activity profile: share of detection minutes in each hour relative to sunrise
    # (bin k = from k h to k+1 h after sunrise; k < 0 = before sunrise)
    dm <- merge(det_min, refs[, .(recorder, date, sunrise_min)], by = c("recorder", "date"))
    dm[, rel_h := floor((m - sunrise_min) / 60)]
    dm[rel_h < -6, rel_h := rel_h + 24L]; dm[rel_h > 17, rel_h := rel_h - 24L]
    bins <- -6:17
    species <- lapply(sort(unique(d_raw$species)), function(sp_name) {
      x <- d_raw[species == sp_name]
      prof <- dm[species == sp_name, .N, by = rel_h][data.table(rel_h = bins), on = "rel_h"]
      prof[is.na(N), N := 0]
      pd <- p_day[species == sp_name]
      bs <- by_sched[species == sp_name]
      list(common_name = sp_name,
           n_detections = nrow(x),
           days_detected_24h = uniqueN(x$date),
           days_detected_by_schedule = setNames(as.list(vapply(targets, function(t)
             if (is.null(rec_rows[[match(t, targets)]])) NA_integer_ else
               as.integer(sum(bs[target == t]$days)), 1L)), paste0(targets, "%")),
           n_detection_minutes = nrow(unique(x[, .(recorder, date, minute_of_day)])),
           confidence = list(mean = round(mean(x$conf), 3), median = round(median(x$conf), 3),
                             max = round(max(x$conf), 3)),
           p_day_detected_by_cycle = setNames(as.list(round(pd$p[match(periods, pd$period)], 3)), paste0(dd, "/", periods)),
           activity_share_by_hour_from_sunrise = round(prof$N / max(1, sum(prof$N)), 3))
    })
  }

  list(
    schema = list(name = "SCORE preset", version = SCORE_PRESET_SCHEMA),
    preset_id = NULL, title = NULL,
    created = format(Sys.Date()), score_version = score_version,
    contributors = NULL, study = NULL,
    sampling = list(
      period = list(first_date = format(min(dates)), last_date = format(max(dates)), n_days = n_days, timing = NULL,
                    mean_sunrise = minutes_to_hm(mean(refs$sunrise_min)), timezone = settings$timezone),
      recorders = list(n_units = settings$n_units, model = NULL, sample_rate_hz = NULL,
                       deployment_height_m = NULL),
      raw_recording = NULL),
    classifier = list(name = NULL, version = NULL, min_confidence = settings$min_conf,
                      species_filter = settings$species_filter),
    analysis = list(spatial_pooling = if (n_rec == 1) "all recorders pooled per day" else "per recorder schedules, pooled richness per day",
                    richness_metric = "mean daily species richness",
                    reference = list(continuous_24h_mean_daily_richness = round(ref24, 3), total_species = uniqueN(dt$species)),
                    duty_cycles = list(record_min = dd, cycle_min = range(sapply(schedules, attr, "period"))),
                    window_grid = settings$window_grid,
                    n_candidates = nrow(combo),
                    selection_rule = "least effort with richness >= target (% of best candidate); ties: highest richness"),
    results = list(recommended_schedules = recs, pareto_front = pareto),
    validation = NULL,
    species = species,
    notes = NULL
  )
}

# Adds the form fields (contributor, study, recording, classifier) to a preset core
# Location: the coordinates entered in the app (Tab 1) are shared only if
# meta$share_location is TRUE, rounded to 0.1 degree (~10 km) if meta$round_location is TRUE.
add_preset_meta <- function(core, meta, settings) {
  p <- core
  location <- NULL
  if (isTRUE(meta$share_location) && !is.null(settings$lat) && !is.na(settings$lat)) {
    prec <- if (isTRUE(meta$round_location)) 0.1 else NULL
    location <- list(lat = if (is.null(prec)) round(settings$lat, 5) else round_coord(settings$lat, prec),
                     lon = if (is.null(prec)) round(settings$lon, 5) else round_coord(settings$lon, prec),
                     rounded = !is.null(prec), precision_deg = prec)
  }
  p$preset_id <- meta$preset_id; p$title <- meta$title
  p$contributors <- list(list(name = meta$name, affiliation = meta$affiliation, orcid = meta$orcid, contact = meta$contact))
  p$study <- list(project_description = meta$project_description, objective = meta$objective,
                  target_taxa = as.list(meta$target_taxa), target_species = as.list(meta$target_species),
                  country = meta$country, country_code = unname(COUNTRIES[meta$country]), region = meta$region,
                  habitat = list(classification = "EUNIS level 1", codes = as.list(meta$habitat_codes), description = meta$habitat_desc),
                  elevation_mean_m = meta$elevation_mean,
                  location = location)
  p$sampling$period$timing <- as.list(meta$timing)
  p$sampling$recorders$model <- meta$recorder_model; p$sampling$recorders$sample_rate_hz <- meta$sample_rate
  p$sampling$recorders$deployment_height_m <- meta$height
  p$sampling$raw_recording <- meta$raw_recording
  p$classifier$name <- meta$classifier; p$classifier$version <- meta$classifier_version
  p$notes <- meta$notes
  # keep a fixed order of the top-level blocks
  p[c("schema", "preset_id", "title", "created", "score_version", "contributors", "study",
      "sampling", "classifier", "analysis", "results", "validation", "species", "notes")]
}

write_score_preset <- function(preset, file) {
  writeLines(jsonlite::toJSON(preset, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA), file)
  invisible(file)
}

# Writes the preset as a .zip: <preset_id>/<preset_id>.json + <preset_id>/figures/*.png
# plots: named list of ggplot objects (figures of the app); NULL elements are skipped
APP_FIGURES <- c(
  duty_cycle_rarefaction = "Tab 2 - Richness by duty cycle (bootstrap)",
  duty_cycle_missed_species = "Tab 2 - Species missed by duty cycle",
  time_window_w1 = "Tab 3 - Richness by time window (window 1)",
  time_window_w2 = "Tab 3 - Richness by time window (window 2)",
  time_window_combined = "Tab 3 - Richness of combined windows",
  multi_optimum = "Tab 4 - Multi-optimum: effort vs richness of all candidate schedules")

write_score_preset_zip <- function(preset, plots, file) {
  id  <- preset$preset_id
  tmp <- tempfile("preset_"); dir.create(file.path(tmp, id, "figures"), recursive = TRUE)
  write_score_preset(preset, file.path(tmp, id, paste0(id, ".json")))
  for (nm in names(plots)) {
    g <- plots[[nm]]; if (is.null(g)) next
    h <- if (nm == "duty_cycle_missed_species")
      min(20, 3 + 0.13 * length(unique(g$data$species_missed))) else 5
    tryCatch(suppressWarnings(ggplot2::ggsave(file.path(tmp, id, "figures", paste0(nm, ".png")), g,
                                              width = 8, height = h, dpi = 120, bg = "white")),
             error = function(e) NULL)
  }
  file <- normalizePath(file, mustWork = FALSE)
  old <- setwd(tmp); on.exit(setwd(old), add = TRUE)
  utils::zip(zipfile = file, files = id, flags = "-r9Xq")
  invisible(file)
}
