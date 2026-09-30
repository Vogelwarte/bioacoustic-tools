# ==============================================================================
# pheno_matrix()  -  phenology heatmap of BirdNET detections
# ------------------------------------------------------------------------------
# Original code: Amandine Serrurier & Jean-Nicolas Pradervand (with C. Sahli).
#
# Counts detections per time slot ('Unit' minutes) and per day, and returns a
# ggplot heatmap with night / twilight in the background. Each tile carries a
# 'text' aesthetic (date, time slot, number of detections) so the plot can be
# made interactive with  plotly::ggplotly(p, tooltip = "text").
#
# IMPORTANT for plotly: the tiles are drawn on the COMPLETE grid (every day x
# every time slot). Empty cells get NA (transparent). If only the non-empty
# cells were drawn, ggplotly would build an irregular heatmap and stretch
# cells across the missing time slots (the long vertical bars).
#
# Arguments
#   Voc          detections, already filtered (needs Start_segment [POSIXct]
#                and Common.Name)
#   SP           species shown in the title ("All species" or a vector)
#   Unit         aggregation interval in minutes
#   Confidence1  kept for compatibility (filtering is done in the app)
#   sunrise      TRUE to draw night / twilight in the background
#   LAT, LONG    coordinates used for the sun times
#   TimeZone     time zone of the time axis (recorder clock or local time)
#   xlim_plot    optional c(start_date, end_date)
#   nocturnal    FALSE = 00:00-24:00 axis, TRUE = 12:00-12:00 axis (one column = one night)
#   tile_alpha   opacity of the tiles (1 = opaque)
#   log_colour   TRUE = log colour scale (rare and frequent slots both readable)
#   fixed_recorder_clock  legacy: shifts sun times by -1 h during DST. Keep
#                FALSE when TimeZone already matches the data.
#
# Returns a ggplot object, or a character string explaining why nothing
# could be plotted.
# ==============================================================================


# --- Colours ------------------------------------------------------------------
# Detections: one-hue orange ramp, light -> dark. Checked: lightness increases
# step by step, adjacent steps are distinguishable, and the lightest step is
# still visible on a white background (the old pale yellow was not).
PHENO_RAMP     <- c("#f0a064", "#e97d3f", "#d45a22", "#a84118", "#72290d", "#401505")
# Background: night = dark blue, dawn / dusk twilight = transparent light blue.
# They are always drawn BELOW the tiles (see pheno_add_sun_shapes()).
PHENO_NIGHT    <- "rgba(27, 42, 94, 0.80)"     # dark blue (plotly)
PHENO_TWILIGHT <- "rgba(120, 180, 235, 0.45)"  # light blue, transparent (plotly)
PHENO_NIGHT_GG    <- "#1B2A5ECC"               # same colours for ggplot (#RRGGBBAA)
PHENO_TWILIGHT_GG <- "#78B4EB73"


pheno_matrix <- function(
    Voc,
    SP = "All species",
    Unit,
    Confidence1 = NULL,
    sunrise = TRUE,
    LAT,
    LONG,
    TimeZone,
    xlim_plot = NULL,
    nocturnal = FALSE,
    tile_alpha = 1,
    log_colour = TRUE,
    fixed_recorder_clock = FALSE,
    sun_as_shapes = FALSE   # TRUE for plotly: sun bands are added later as shapes
) {
  
  # --- 1. INPUT DATA -----------------------------------------------------------
  dt <- data.table::as.data.table(data.table::copy(Voc))
  if (nrow(dt) == 0) {
    return("No detections for this selection (species, recorders or confidence threshold).")
  }
  
  # Timestamps expressed in the plotting time zone
  dt[, Start_local := lubridate::with_tz(Start_segment, tzone = TimeZone)]
  
  # Plot title
  n_sp <- length(unique(dt$Common.Name))
  SP_title <- if (length(SP) > 1) {
    paste0("Several species selected (n = ", n_sp, ")")
  } else if (SP == "All species") {
    paste0("All detected species (n = ", n_sp, ")")
  } else {
    SP
  }
  
  # --- 2. COUNTS PER DAY AND TIME SLOT ------------------------------------------
  slot_h  <- Unit / 60                                  # slot height in hours
  n_slots <- ceiling(1440 / Unit)
  
  dt[, date := as.Date(Start_local, tz = TimeZone)]
  dt[, slot := pmin(floor((lubridate::hour(Start_local) * 60 +
                             lubridate::minute(Start_local)) / Unit), n_slots - 1)]
  counts <- dt[, .N, by = .(date, slot)]
  counts[, slot_start := slot * slot_h]                 # decimal hours, 0 .. 24
  
  # Nocturnal axis: morning slots (< 12 h) belong to the previous night, +24 h
  if (nocturnal) {
    counts[, x := data.table::fifelse(slot_start < 12, date - 1, date)]
    counts[, y := data.table::fifelse(slot_start < 12, slot_start + 24, slot_start)]
    y_lim <- c(12, 36)
  } else {
    counts[, x := date]
    counts[, y := slot_start]
    y_lim <- c(0, 24)
  }
  
  # --- 3. DATE RANGE ------------------------------------------------------------
  x_range <- if (!is.null(xlim_plot)) as.Date(xlim_plot) else range(counts$x)
  counts  <- counts[x >= x_range[1] & x <= x_range[2]]
  if (nrow(counts) == 0) return("No data to plot in the selected date range.")
  
  # --- 4. COMPLETE GRID (every day x every slot) --------------------------------
  grid <- data.table::CJ(x = seq(x_range[1], x_range[2], by = "day"),
                         y = round(seq(y_lim[1], y_lim[2] - slot_h, by = slot_h), 6))
  counts[, y := round(y, 6)]
  tiles <- merge(grid, counts[, .(x, y, N)], by = c("x", "y"), all.x = TRUE)
  # N = NA where there is no detection -> transparent tile, no tooltip
  
  # Tooltip: real calendar date and time slot of each cell
  fmt_hm    <- function(m) sprintf("%02d:%02d", (m %/% 60) %% 24, m %% 60)
  start_min <- round((tiles$y %% 24) * 60)
  real_date <- tiles$x + as.integer(tiles$y >= 24)       # nocturnal: after midnight = next day
  tiles[, text := data.table::fifelse(
    is.na(N), "",
    paste0("Date: ", format(real_date, "%a %d %b %Y"),
           "<br>Time: ", fmt_hm(start_min), " - ", fmt_hm(start_min + Unit),
           "<br>Detections: ", N)
  )]
  tiles[, y_mid := y + slot_h / 2]                        # centre of the slot
  
  # --- 5. SUN TIMES -------------------------------------------------------------
  SunDF <- NULL
  if (sunrise) {
    to_hour <- function(x) lubridate::hour(x) + lubridate::minute(x) / 60
    get_sun <- function(d) {
      S <- suncalc::getSunlightTimes(date = d, lat = LAT, lon = LONG,
                                     keep = c("sunrise", "sunset", "dawn", "dusk"),
                                     tz = TimeZone)
      if (fixed_recorder_clock) {
        shift <- lubridate::hours(ifelse(lubridate::dst(S$sunrise), 1, 0))
        for (v in c("sunrise", "sunset", "dawn", "dusk")) S[[v]] <- S[[v]] - shift
      }
      S
    }
    days <- seq(x_range[1], x_range[2], by = "day")
    SunDF <- tryCatch({
      if (!nocturnal) {
        S <- get_sun(days)
        data.frame(x = days, dawn = to_hour(S$dawn), sunrise = to_hour(S$sunrise),
                   sunset = to_hour(S$sunset), dusk = to_hour(S$dusk))
      } else {
        S_eve <- get_sun(days)       # evening of day D
        S_mor <- get_sun(days + 1)   # morning of day D + 1 (same night)
        data.frame(x = days, sunset = to_hour(S_eve$sunset), dusk = to_hour(S_eve$dusk),
                   dawn = to_hour(S_mor$dawn) + 24, sunrise = to_hour(S_mor$sunrise) + 24)
      }
    }, error = function(e) {
      warning("Sun time calculation failed: ", e$message)
      NULL
    })
    if (!is.null(SunDF)) SunDF <- SunDF[stats::complete.cases(SunDF), ]
    if (!is.null(SunDF) && nrow(SunDF) == 0) SunDF <- NULL
  }
  
  # --- 6. PLOT ------------------------------------------------------------------
  a <- ggplot2::ggplot()
  
  # Background: night and twilight.
  # For plotly (sun_as_shapes = TRUE) they are NOT drawn here: plotly always
  # draws ribbons above heatmaps, so the app adds them as shapes below the tiles.
  if (!is.null(SunDF) && !sun_as_shapes) {
    SunDF$bottom <- y_lim[1]
    SunDF$top    <- y_lim[2]
    band <- function(lo, hi, col) {
      ggplot2::geom_ribbon(data = SunDF,
                           ggplot2::aes(x = x, ymin = .data[[lo]], ymax = .data[[hi]]),
                           fill = col, colour = NA)
    }
    if (!nocturnal) {
      a <- a + band("bottom", "dawn", PHENO_NIGHT_GG) + band("dawn", "sunrise", PHENO_TWILIGHT_GG) +
        band("sunset", "dusk", PHENO_TWILIGHT_GG) + band("dusk", "top", PHENO_NIGHT_GG)
    } else {
      a <- a + band("sunset", "dusk", PHENO_TWILIGHT_GG) + band("dusk", "dawn", PHENO_NIGHT_GG) +
        band("dawn", "sunrise", PHENO_TWILIGHT_GG)
    }
  }
  
  # Tiles on the complete grid (NA = transparent)
  use_log <- log_colour && max(tiles$N, na.rm = TRUE) > 1
  a <- a +
    ggplot2::geom_tile(
      data = tiles,
      ggplot2::aes(x = x, y = y_mid, fill = N, text = text),
      width = 1, height = slot_h, alpha = tile_alpha
    ) +
    ggplot2::scale_fill_gradientn(
      colours  = PHENO_RAMP,
      trans    = if (use_log) "log10" else "identity",
      limits   = c(1, max(2, max(tiles$N, na.rm = TRUE))),
      na.value = "transparent",
      name     = paste0("Detections /\n", Unit, " min")
    ) +
    ggplot2::scale_y_continuous(
      "Time",
      limits = y_lim,
      breaks = seq(y_lim[1], y_lim[2], by = 2),
      labels = function(v) sprintf("%02d:00", as.integer(v %% 24)),
      expand = c(0, 0)
    ) +
    ggplot2::scale_x_date("Date", date_labels = "%d %b", expand = c(0, 0)) +
    ggplot2::labs(title = SP_title) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white"),
      plot.background  = ggplot2::element_rect(fill = "white"),
      panel.grid       = ggplot2::element_blank(),
      text             = ggplot2::element_text(size = 14)
    )
  
  # Keep the sun times with the plot so the app can draw them as plotly shapes
  attr(a, "sun")       <- SunDF
  attr(a, "nocturnal") <- nocturnal
  attr(a, "y_lim")     <- y_lim
  a
}


# ==============================================================================
# pheno_add_sun_shapes()  -  night / twilight bands BELOW the tiles (plotly)
# ------------------------------------------------------------------------------
# plotly always draws scatter / ribbon traces above heatmap traces, whatever
# their order. Layout shapes with layer = "below" are drawn under all traces,
# so the detections stay on top. One rectangle per day and per band.
#
#   pl  plotly object from ggplotly(p)
#   p   the ggplot returned by pheno_matrix(..., sun_as_shapes = TRUE)
# ==============================================================================
pheno_add_sun_shapes <- function(pl, p) {
  for (i in seq_along(pl$x$data)) {
    tr <- pl$x$data[[i]]
    if (identical(tr$type, "heatmap")) {
      cs <- tr$colorscale
      if (!is.null(cs) && NROW(cs) < 2) {
        col <- if (is.data.frame(cs)) cs[1, 2] else cs[[1]][[2]]
        pl$x$data[[i]]$colorscale <- list(c(0, col), c(1, col))
      }
    }
  }
  sun <- attr(p, "sun")
  if (is.null(sun) || nrow(sun) == 0) return(pl)
  noct  <- isTRUE(attr(p, "nocturnal"))
  y_lim <- attr(p, "y_lim")
  
  # x coordinates: ggplotly may use a date axis or numeric days - handle both
  x_is_date <- identical(pl$x$layout$xaxis$type, "date")
  x_left  <- sun$x - 0.5
  x_right <- sun$x + 0.5
  to_x <- function(d) {
    if (x_is_date) format(as.POSIXct(as.numeric(d) * 86400, origin = "1970-01-01", tz = "UTC"),
                          "%Y-%m-%d %H:%M:%S")
    else as.numeric(d)
  }
  x0 <- to_x(x_left); x1 <- to_x(x_right)
  
  rect <- function(i, ylo, yhi, col) {
    list(type = "rect", layer = "below", xref = "x", yref = "y",
         x0 = x0[i], x1 = x1[i], y0 = ylo, y1 = yhi,
         fillcolor = col, line = list(width = 0))
  }
  
  shapes <- list()
  for (i in seq_len(nrow(sun))) {
    s <- sun[i, ]
    if (!noct) {
      shapes <- c(shapes, list(
        rect(i, y_lim[1],  s$dawn,    PHENO_NIGHT),
        rect(i, s$dawn,    s$sunrise, PHENO_TWILIGHT),
        rect(i, s$sunset,  s$dusk,    PHENO_TWILIGHT),
        rect(i, s$dusk,    y_lim[2],  PHENO_NIGHT)
      ))
    } else {
      shapes <- c(shapes, list(
        rect(i, s$sunset,  s$dusk,    PHENO_TWILIGHT),
        rect(i, s$dusk,    s$dawn,    PHENO_NIGHT),
        rect(i, s$dawn,    s$sunrise, PHENO_TWILIGHT)
      ))
    }
  }
  
  pl$x$layout$shapes <- c(pl$x$layout$shapes, shapes)
  pl$x$layout$plot_bgcolor <- "white"
  pl
}
