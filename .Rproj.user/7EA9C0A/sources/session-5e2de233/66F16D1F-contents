# ---- Helper Utilities ----

# ---- Input validation ----


# Check that all required column names are present in a data.table.
# Raises a Shiny validation error listing any missing columns.
need_cols <- function(dt, cols, context = "data") {
  missing <- setdiff(cols, names(dt))
  validate(
    need(
      length(missing) == 0,
      paste0("Missing required columns in ", context, ": ", paste(missing, collapse = ", "))
    )
  )
  invisible(TRUE)
}

# ---- File path utilities ----

# Normalise BirdNET Begin.Path: convert backslashes to forward slashes and
# collapse runs of slashes to a single slash.
normalize_begin_path <- function(x) {
  x <- gsub("\\\\+", "/", x)
  x <- gsub("/+",    "/", x)
  x
}

# Read a recorder coordinate file.
# Supports .csv, .txt (tab-separated with header), and .xlsx.
# Column names are lowercased before returning.
# Returns NULL for unsupported extensions.
read_coord_file <- function(path, name) {
  ext <- tools::file_ext(name)
  if (ext == "csv") {
    df <- read.csv(path)
  } else if (ext == "txt") {
    df <- read.table(path, header = TRUE)
  } else if (ext == "xlsx") {
    df <- readxl::read_excel(path)
  } else {
    return(NULL)
  }
  names(df) <- tolower(names(df))
  df
}

# Split a Begin.Path column into a matrix of path components.
# Each column Vi contains the i-th segment of the path.
# Returns NULL if x is empty.
split_path_matrix <- function(x, sep = "/") {
  if (length(x) == 0) return(NULL)
  parts   <- stringr::str_split(x, sep)
  max_len <- max(lengths(parts))
  mat <- t(vapply(parts, function(p) {
    length(p) <- max_len
    p
  }, FUN.VALUE = character(max_len)))
  colnames(mat) <- paste0("V", seq_len(ncol(mat)))
  mat
}

# Lightweight preview version of split_path_matrix: only parses the first
# n_preview rows, to avoid heavy computation during UI rendering.
split_path_matrix_preview <- function(x, sep = "/", n_preview = 20L) {
  if (length(x) == 0) return(NULL)
  x <- x[seq_len(min(length(x), n_preview))]
  split_path_matrix(x, sep = sep)
}

#Erase filters

# ---- Time-window utilities ----

# Split a recording window that may span midnight into at most two segments,
# each anchored to the correct calendar date.
#
# Arguments:
#   date      – reference Date
#   start_min – window start in minutes from midnight (may exceed 0–1440)
#   end_min   – window end   in minutes from midnight (may exceed 0–1440)
#
# Returns a data.table with columns: date, start_min, end_min.
split_window_across_days <- function(date, start_min, end_min) {
  
  s_shift <- start_min %/% 1440
  e_shift <- end_min   %/% 1440
  
  s_mod <- start_min %% 1440
  e_mod <- end_min   %% 1440
  
  # Both endpoints fall on the same relative day
  if (s_shift == e_shift) {
    return(data.table(
      date      = as.Date(date) + s_shift,
      start_min = s_mod,
      end_min   = e_mod
    ))
  }
  
  # Window crosses midnight: two segments
  data.table(
    date      = c(as.Date(date) + s_shift, as.Date(date) + e_shift),
    start_min = c(s_mod, 0L),
    end_min   = c(1440L, e_mod)
  )
}

# Test whether a recording block [start_min, end_min] falls within a recording
# window [win_start, win_end]. Handles windows that overflow into the next day
# (win_end > 1440) or the previous day (win_start < 0).
#
# win_start and win_end are treated as scalars (first element used).
# start_min / end_min may be vectors.
block_in_window <- function(start_min, end_min, win_start, win_end) {
  
  # Enforce scalar window bounds
  win_start <- win_start[1]
  win_end   <- win_end[1]
  
  # Normal case: window stays within a single day
  if (win_start >= 0 && win_end <= 1440) {
    return(start_min >= win_start & end_min <= win_end)
  }
  
  # Window overflows into the next day
  if (win_end > 1440) {
    return(
      (start_min >= pmax(0, win_start) & end_min <= 1440) |
        (start_min >= 0 & end_min <= (win_end - 1440))
    )
  }
  
  # Window starts before midnight of the previous day
  if (win_start < 0) {
    return(
      (start_min >= (1440 + win_start) & end_min <= 1440) |
        (start_min >= 0 & end_min <= pmin(1440, win_end))
    )
  }
  
  return(rep(FALSE, length(start_min)))
}



# ---- Schedule generation ----

# Build a duty-cycle recording schedule aligned to solar events.
#
# For each date in sun_table the function computes an active recording window
# (defined by a solar start/end event with optional hour extensions or a fixed
# duration) and then tiles that window with duty-cycle blocks.
#
# Arguments:
#   sun_table           – data.table with columns: date, dawn, sunrise, sunset, dusk
#   start_event         – column name of the window start solar event
#   end_event           – column name of the window end   solar event
#   extend_before_hours – hours to shift the window start earlier
#   extend_after_hours  – hours to shift the window end   later
#   length_morning      – if > 0, use a fixed window duration (hours) instead of end_event
#   period              – duty-cycle period in minutes
#   duty_duration       – recording block length in minutes (must be ≤ period)
#   nb_duty             – number of duty blocks per period (currently kept for API
#                         compatibility; expansion handled by the caller)
#   random_start        – if TRUE, randomise the block start within each period
#   record_24h          – if TRUE, record continuously (00:00–24:00) regardless of events
#
# Returns a data.table with columns: date, start_min, end_min.
create_sun_duty_schedule <- function(sun_table,
                                     start_event         = "dawn",
                                     end_event           = "dusk",
                                     extend_before_hours = 1,
                                     extend_after_hours  = 3,
                                     length_morning      = 0,
                                     period              = 5,
                                     duty_duration       = 1,
                                     nb_duty             = 1,
                                     random_start        = FALSE,
                                     record_24h          = FALSE) {
  
  # -- 1. Validation --
  period        <- as.numeric(period)
  duty_duration <- as.numeric(duty_duration)
  
  if (duty_duration > period) {
    stop("duty_duration cannot be greater than period.")
  }
  
  DT <- as.data.table(sun_table)
  
  # -- 2. Compute recording window boundaries (vectorised over all dates) --
  if (record_24h) {
    DT[, `:=`(
      start_min = 0L,
      end_min   = 1440L,
      ref_date  = as.Date(get(start_event))
    )]
  } else {
    start_times <- get(start_event, DT) - dhours(extend_before_hours)
    
    if (!is.null(length_morning) && length_morning > 0) {
      end_times <- start_times + dhours(length_morning)
    } else {
      end_times <- get(end_event, DT) + dhours(extend_after_hours)
    }
    
    DT[, `:=`(
      start_min = as.integer(hour(start_times) * 60 + minute(start_times)),
      end_min   = as.integer(hour(end_times)   * 60 + minute(end_times)),
      ref_date  = as.Date(start_times)
    )]
  }
  
  # -- 3. Generate duty-cycle block start times for each date --
  generate_blocks <- function(s_min, e_min, date_ref) {
    if (s_min > e_min) {
      # Window crosses midnight: two segments
      part1 <- seq(s_min, 1439, by = period)
      part2 <- seq(0, e_min, by = period)
      blocks <- c(part1, part2)
      dates  <- c(rep(date_ref, length(part1)), rep(date_ref + 1, length(part2)))
    } else {
      blocks <- seq(s_min, e_min, by = period)
      dates  <- rep(date_ref, length(blocks))
    }
    
    if (length(blocks) == 0) return(NULL)
    
    data.table(
      date         = as.Date(dates),
      block_start  = as.integer(blocks),
      window_start = s_min,
      window_end   = e_min
    )
  }
  
  blocks_list <- mapply(
    generate_blocks,
    DT$start_min,
    DT$end_min,
    DT$ref_date,
    SIMPLIFY = FALSE
  )
  
  all_blocks <- rbindlist(blocks_list, fill = TRUE)
  if (nrow(all_blocks) == 0) return(data.table())
  
  # -- 4. Expand for nb_duty > 1 and compute start_min / end_min --
  if (nb_duty > 1) {
    all_blocks <- all_blocks[rep(seq_len(.N), each = nb_duty)]
  }
  
  if (random_start) {
    max_offset <- period - duty_duration
    offsets <- mapply(function(limit) {
      if (limit < 0) return(0)
      sample(0:limit, 1)
    }, max_offset)
    all_blocks[, start_min := block_start + offsets]
  } else {
    all_blocks[, start_min := block_start]
  }
  
  all_blocks[, end_min := start_min + duty_duration - 1L]
  
  # Remove intermediate columns
  all_blocks[, `:=`(block_start = NULL, window_start = NULL, window_end = NULL)]
  
  return(all_blocks[])
}


# ---- Species richness calculation ----

# Compute species richness from a detection table filtered by a duty-cycle
# schedule.
#
# Arguments:
#   dt           – detection data.table with columns: date, minute_of_day, species
#                  (optionally: recorder, site)
#   schedule     – schedule data.table with columns: date, start_min, end_min
#                  (optionally: recorder)
#   per_spatial  – spatial grouping: "recorder", "site", or "all"
#   per_temporal – temporal grouping: "day", "month", or "season"
#   return_species – if TRUE, also return the list of detected species
#
# Returns a data.frame of richness (and effort in minutes) by group,
# or a list(richness, species) when return_species = TRUE.
compute_richness <- function(dt, schedule,
                             per_spatial  = c("recorder", "site", "all"),
                             per_temporal = c("day", "month", "season"),
                             return_species = FALSE) {
  
  per_spatial  <- match.arg(per_spatial)
  per_temporal <- match.arg(per_temporal)
  
  # -- 1. Type safety --
  if (!is.data.table(dt))       dt       <- as.data.table(dt)
  if (!is.data.table(schedule)) schedule <- as.data.table(schedule)
  
  if (nrow(schedule) == 0) stop("Schedule is empty.")
  if (!"minute_of_day" %in% names(dt))       stop("Column 'minute_of_day' missing in dt.")
  if (!"date"          %in% names(dt))       stop("Column 'date' missing in dt.")
  if (!"date"          %in% names(schedule)) stop("Column 'date' missing in schedule.")
  
  # Force Date type consistently
  dt[,       date := as.Date(date, origin = "1970-01-01")]
  schedule[, date := as.Date(date, origin = "1970-01-01")]
  
  dt[,       minute_of_day := as.integer(minute_of_day)]
  schedule[, `:=`(start_min = as.integer(start_min), end_min = as.integer(end_min))]
  
  # -- 2. Resolve spatial column conflicts between dt and schedule --
  schedule_work <- copy(schedule)
  
  has_recorder_dt    <- "recorder" %in% names(dt)
  has_recorder_sched <- "recorder" %in% names(schedule_work)
  
  if (has_recorder_sched && !has_recorder_dt) {
    setnames(schedule_work, "recorder", "recorder_ignore")
  }
  
  # -- 3. Non-equi join: keep only detections that fall inside a scheduled block --
  join_cols <- "date"
  if (has_recorder_dt && "recorder" %in% names(schedule_work)) {
    join_cols <- c(join_cols, "recorder")
  }
  
  on_conditions <- c(
    paste0(join_cols, "==", join_cols),
    "start_min<=minute_of_day",
    "end_min>=minute_of_day"
  )
  
  dt_filt <- schedule_work[dt, on = on_conditions, nomatch = 0L]
  
  # -- 4. Keep only original dt columns --
  cols_to_keep <- names(dt)
  existing_keep <- cols_to_keep[cols_to_keep %in% names(dt_filt)]
  
  if (length(existing_keep) == 0) {
    tech_cols    <- c("start_min", "end_min", "recorder_ignore")
    existing_keep <- names(dt_filt)[!names(dt_filt) %in% tech_cols]
  }
  
  dt_filt <- dt_filt[, ..existing_keep]
  
  # -- 5. Compute recording effort per day (from the original schedule) --
  if ("recorder" %in% names(schedule)) {
    effort_daily <- schedule[, .(effort_minutes = sum(end_min - start_min + 1)),
                             by = .(recorder, date)]
  } else {
    effort_daily <- schedule[, .(effort_minutes = sum(end_min - start_min + 1)),
                             by = .(date)]
  }
  effort_daily[, date := as.Date(date, origin = "1970-01-01")]
  
  if (nrow(dt_filt) == 0) {
    warning("No detections found within the scheduled windows.")
    return(NULL)
  }
  
  df_work <- as.data.frame(dt_filt)
  
  # -- 6. Add temporal grouping columns --
  if (per_temporal == "month") {
    df_work$month <- format(df_work$date, "%Y-%m")
  } else if (per_temporal == "season") {
    mois  <- month(df_work$date)
    annee <- year(df_work$date)
    df_work$season_year <- ifelse(
      mois %in% c(12, 1, 2),
      paste0("Winter_",  ifelse(mois == 12, annee, annee - 1), "-", ifelse(mois == 12, annee + 1, annee)),
      ifelse(mois %in% c(3, 4, 5),  paste0("Spring_", annee),
             ifelse(mois %in% c(6, 7, 8), paste0("Summer_", annee),
                    paste0("Autumn_", annee)))
    )
  }
  
  # -- 7. Build grouping vector --
  spatial_col <- NULL
  if (per_spatial == "recorder" && "recorder" %in% names(df_work)) spatial_col <- "recorder"
  if (per_spatial == "site"     && "site"     %in% names(df_work)) spatial_col <- "site"
  
  temporal_col <- switch(per_temporal,
                         day    = "date",
                         month  = "month",
                         season = "season_year")
  
  groupes <- c(spatial_col, temporal_col)
  groupes <- groupes[!sapply(groupes, is.null)]
  
  if (length(groupes) == 0) {
    df_work$groupe_global <- "Total"
    groupes <- "groupe_global"
  }
  
  # -- 8. Compute richness --
  resultat <- df_work %>%
    group_by(across(all_of(groupes))) %>%
    summarise(richness = n_distinct(species), .groups = "drop")
  
  # -- 9. Append effort --
  if (!is.null(effort_daily) && length(groupes) > 0) {
    effort_df <- as.data.frame(effort_daily)
    
    if ("month"       %in% groupes) effort_df$month       <- format(effort_df$date, "%Y-%m")
    if ("season_year" %in% groupes) {
      mois  <- month(effort_df$date)
      annee <- year(effort_df$date)
      effort_df$season_year <- ifelse(
        mois %in% c(12, 1, 2),
        paste0("Winter_",  ifelse(mois == 12, annee, annee - 1), "-", ifelse(mois == 12, annee + 1, annee)),
        ifelse(mois %in% c(3, 4, 5),  paste0("Spring_", annee),
               ifelse(mois %in% c(6, 7, 8), paste0("Summer_", annee),
                      paste0("Autumn_", annee)))
      )
    }
    if ("groupe_global" %in% groupes) effort_df$groupe_global <- "Total"
    
    effort_agg <- effort_df %>%
      group_by(across(all_of(groupes))) %>%
      summarise(effort_minutes = sum(effort_minutes), .groups = "drop")
    
    resultat <- left_join(resultat, effort_agg, by = groupes)
  }
  
  if ("groupe_global" %in% names(resultat)) resultat$groupe_global <- NULL
  
  # -- 10. Optionally return species list alongside richness --
  if (return_species) {
    liste_especes <- df_work %>%
      group_by(across(all_of(groupes)), species) %>%
      summarise(.groups = "drop")
    return(list(richness = resultat, species = liste_especes))
  }
  
  return(resultat)
}

