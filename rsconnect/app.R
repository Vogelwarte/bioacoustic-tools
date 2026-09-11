
#Amandine SERRURIER 13.01.2026
# function to rarify recording and accumulate richness AND corresponding shiny app

# ### TOKEN RSCONNECT#####
#  setwd('/Users/ase/Library/Mobile Documents/com~apple~CloudDocs/Documents/VoWa/shiny_App_raref')
# # # # #
# library(rsconnect)
# rsconnect::forgetDeployment()
# rsconnect::setAccountInfo(name='vogelwarte',
#                           token='5FA01481ECE5F33483C37F69731D3CF9',
#                           secret='bBVKSYJlP/1Xz2Ynz8FuwxotsbgoQQTRrKekRXh7')

# deployApp()


# ---- Libraries ----
library(shiny)
library(shinyFiles)
library(shinythemes)
library(data.table)
library(lubridate)
library(ggplot2)
library(future.apply)
library(suncalc)
library(tidyr)
library(dplyr)
library(DT)
library(scico)
library(plotly)
library(colorspace)
library(stringr)

# ---- Global Options ----
options(future.globals.maxSize = 8 * 1024^3)  # 8 GB memory limit for future processing

# ---- Helper Utilities ----

# Lightweight debug logger (console only)
log_debug <- function(...) {
  if (isTRUE(getOption("app.debug"))) message(Sys.time(), ...)
}

# Validate required columns early with consistent UX
need_cols <- function(dt, cols, context = "data") {
  missing <- setdiff(cols, names(dt))
  validate(
    need(length(missing) == 0,
         paste0("Missing required columns in ", context, ": ", paste(missing, collapse = ", ")))
  )
  invisible(TRUE)
}

# Normalize BirdNET Begin.Path safely (same transformations as before)
normalize_begin_path <- function(x) {
  x <- gsub("\\\\+", "/", x)
  x <- gsub("/+", "/", x)
  x
}

# Robust coordinate file reader (same behavior; safer errors)
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

validate_windows_against_domain <- function(
    design,
    solar_row,
    domain_start_event,
    domain_end_event,
    grid_event1,
    duration1_col,
    start1_col,
    grid_event2 = NULL,
    duration2_col = NULL,
    start2_col = NULL,
    timezone
) {
  
  # --- helper : convertir event en minute absolue ---
  get_event_min <- function(event_name, date, solar_row, timezone) {
    
    if (event_name %in% c("dawn", "sunrise", "sunset", "dusk")) {
      ref <- solar_row[[event_name]]
    } else if (grepl("^\\d{2}:\\d{2}$", event_name)) {
      ref <- as.POSIXct(paste(date, event_name), tz = timezone)
    } else {
      stop("Unknown event name")
    }
    
    return(lubridate::hour(ref) * 60 + lubridate::minute(ref))
  }
  
  # --- Domaine expérimental absolu ---
  domain_start_min <- get_event_min(domain_start_event,
                                    solar_row$date,
                                    solar_row,
                                    timezone)
  
  domain_end_min   <- get_event_min(domain_end_event,
                                    solar_row$date,
                                    solar_row,
                                    timezone)
  
  # --- Vérification Window 1 ---
  ref1_min <- get_event_min(grid_event1,
                            solar_row$date,
                            solar_row,
                            timezone)
  
  start1_abs <- ref1_min + design[[start1_col]] * 60
  end1_abs   <- start1_abs + design[[duration1_col]] * 60
  
  valid1 <- start1_abs >= domain_start_min &
    end1_abs   <= domain_end_min
  
  # --- Vérification Window 2 (si présente) ---
  if (!is.null(grid_event2)) {
    
    ref2_min <- get_event_min(grid_event2,
                              solar_row$date,
                              solar_row,
                              timezone)
    
    start2_abs <- ref2_min + design[[start2_col]] * 60
    end2_abs   <- start2_abs + design[[duration2_col]] * 60
    
    valid2 <- start2_abs >= domain_start_min &
      end2_abs   <= domain_end_min
    
    return(valid1 & valid2)
    
  } else {
    return(valid1)
  }
}

allowed_offset_range_h <- function(solar_row,
                                   domain_start_event,
                                   domain_end_event,
                                   ref_event,
                                   timezone) {
  
  get_event_min <- function(event_name, date, solar_row, timezone) {
    
    if (event_name %in% c("dawn", "sunrise", "sunset", "dusk")) {
      ref <- solar_row[[event_name]]
    } else if (grepl("^\\d{2}:\\d{2}$", event_name)) {
      ref <- as.POSIXct(paste(date, event_name), tz = timezone)
    } else {
      stop("Unknown event name: ", event_name)
    }
    
    lubridate::hour(ref) * 60 + lubridate::minute(ref)
  }
  
  d <- as.Date(solar_row$date)
  
  domain_start_abs <- get_event_min(domain_start_event, d, solar_row, timezone)
  domain_end_abs   <- get_event_min(domain_end_event,   d, solar_row, timezone)
  ref_abs          <- get_event_min(ref_event,          d, solar_row, timezone)
  
  # offsets relatifs à l'event de référence
  min_h <- (domain_start_abs - ref_abs) / 60
  max_h <- (domain_end_abs   - ref_abs) / 60
  
  c(min_h = min_h, max_h = max_h)
}

split_window_across_days <- function(date, start_min, end_min) {
  # date: Date (jour de référence)
  # start_min / end_min: minutes relatives à ce jour (peuvent dépasser 0..1440)
  
  s_shift <- start_min %/% 1440
  e_shift <- end_min   %/% 1440
  
  s_mod <- start_min %% 1440
  e_mod <- end_min   %% 1440
  
  # même jour "relatif"
  if (s_shift == e_shift) {
    return(data.table(
      date = as.Date(date) + s_shift,
      start_min = s_mod,
      end_min = e_mod
    ))
  }
  
  # traverse minuit -> 2 segments
  # segment 1: de start -> 1440 sur le jour start
  # segment 2: de 0 -> end sur le jour end
  data.table(
    date = c(as.Date(date) + s_shift, as.Date(date) + e_shift),
    start_min = c(s_mod, 0L),
    end_min   = c(1440L, e_mod)
  )
}
# ---- Global Functions ----

# Function: create schedule based on solar events and duty cycle parameters
create_sun_duty_schedule <- function(sun_table,
                                     start_event = "dawn",
                                     end_event   = "dusk",
                                     extend_before_hours = 1,
                                     extend_after_hours  = 3,
                                     length_morning = 0,
                                     period = 5,
                                     duty_duration = 1,
                                     nb_duty = 1,
                                     random_start = FALSE,
                                     record_24h = FALSE) {
  
  # Ensure numeric inputs (original behavior)
  period <- as.numeric(period)
  duty_duration <- as.numeric(duty_duration)
  
  # Step 1: Get sunlight times (sun_table contains: recorder, date, dawn, sunrise, sunset, dusk)
  sun_times <- sun_table
  if (record_24h) {
    
    daily <- tibble(
      start_time = sun_times[[start_event]],
      end_time   = sun_times[[start_event]]
    )
    
    daily$start_time_min <- 0L
    daily$end_time_min   <- 1440L
    
  } else {
    
    # ton code existant
    
    # Determine start and end time based on input conditions
    start_time <- sun_times[[start_event]] - as.period(dhours(extend_before_hours))
    if (!is.null(length_morning) && length_morning > 0) {
      end_time <- start_time + as.period(dhours(length_morning))
    } else {
      end_time <- sun_times[[end_event]] + as.period(dhours(extend_after_hours))
    }
    
    # Step 2: Convert times to minutes of day
    daily <- tibble(start_time, end_time)
    daily$start_time_min <- hour(start_time) * 60 + minute(start_time)
    daily$end_time_min   <- hour(end_time) * 60 + minute(end_time)
    daily$start_time_min <- as.integer(daily$start_time_min)
    daily$end_time_min   <- as.integer(daily$end_time_min)
  }
  # Step 3: Create duty blocks within the defined window
  daily_list <- list()
  window <- list() # kept (original unused var)
  for (i in seq_len(nrow(daily))) {
    day_row <- daily[i, ]
    if (day_row$start_time_min <= day_row$end_time_min) {
      
      blocks <- seq(
        day_row$start_time_min,
        day_row$end_time_min,
        by = period
      )
      
    } else {
      
      # ---- Crossing midnight ----
      
      part1 <- seq(
        day_row$start_time_min,
        1440 - 1,
        by = period
      )
      
      part2 <- seq(
        0,
        day_row$end_time_min,
        by = period
      )
      
      blocks <- c(part1, part2)
    }
    # blocks <- seq(day_row$start_time_min, day_row$end_time_min, by = period)
    
    duty_list <- list()
    daily_window <- list() # kept (original unused var)
    idx <- 1L
    for (b in blocks) {
      for (j in seq_len(nb_duty)) {
        if (random_start) {
          start_duty <- sample(b:(b + period - duty_duration), 1)
        } else {
          start_duty <- b
        }
        end_duty <- start_duty + duty_duration - 1
        duty_list[[idx]] <- data.frame(
          date = as.Date(day_row$start_time),
          start_min = start_duty,
          end_min   = end_duty,
          start_time_min = day_row$start_time_min,
          end_time_min   = day_row$end_time_min
        )
        idx <- idx + 1L
      }
    }
    daily_list[[i]] <- rbindlist(duty_list)
  }
  
  schedule <- do.call(rbind, daily_list)
  schedule$date <- as.Date(schedule$date)
  
  return(schedule)
}

# Function: bootstrap recording units to assess richness accumulation
bootstrap_recorders <- function(dt, B = 100) {
  recorders <- unique(dt$recorder)
  N <- length(recorders)
  k_vals <- seq_len(N)
  
  res <- lapply(k_vals, function(k) {
    lapply(seq_len(B), function(b) {
      sampled_recorders <- sample(recorders, size = k, replace = TRUE)
      data.table(
        k_recorders = k,
        bootstrap   = b,
        recorder    = sampled_recorders
      )
    }) %>% rbindlist()
  }) %>% rbindlist()
  
  return(res)
}

# Function: compute richness per schedule and grouping scheme
# NOTE: Signature preserved; SAFE internal organization only.
compute_richness <- function(dt, schedule,
                             per_spatial = c("recorder", "site", "all"),
                             per_temporal = c("day", "month", "season"),
                             return_species = FALSE) {
  
  per_spatial  <- match.arg(per_spatial)
  per_temporal <- match.arg(per_temporal)
  
  dt$date <- as.Date(dt$date)
  
  if (!is.null(schedule)) {
    
    # Add start/end minute for matching (original behavior: adds columns by reference)
    dt[, `:=`(start_min = minute_of_day, end_min = minute_of_day)]
    setDT(dt)
    setDT(schedule)
    
    # Key (original behavior)
    setkey(schedule, date, start_min, end_min)
    
    # 1. Filter by daily window (non-equi join)
    dt_filt <- schedule[
      dt,
      on = .(recorder, date, start_min <= minute_of_day, end_min >= minute_of_day),
      nomatch = 0
    ]
  } else {
    # spatial case no schedule
    dt_filt <- copy(dt)
  }
  
  # Enrich with month info
  dt_filt <- dt_filt %>% mutate(month = format(date, "%Y-%m"))
  
  if (!is.null(schedule) && nrow(schedule) == 0) {
    stop("Schedule is empty.")
  }
  
  # Calculate daily effort
  if (!is.null(schedule)) {
    effort_by_day <- schedule %>%
      group_by(date) %>%
      summarise(
        effort_minutes = sum(end_min - start_min + 1),
        .groups = "drop"
      )
  } else {
    effort_by_day <- NULL
  }
  
  # Determine spatial grouping
  spatial_group <- switch(per_spatial,
                          recorder = c("recorder"),
                          site = c("site"),
                          all = NULL)
  
  # Determine temporal grouping
  temporal_group <- switch(per_temporal,
                           day = c("date"),
                           month = c("month"),
                           season = NULL)
  
  # Combine grouping levels
  groups <- c(spatial_group, temporal_group)
  
  if (per_temporal %in% c("day", "month")) {
    out <- dt_filt %>%
      group_by(across(all_of(groups))) %>%
      summarise(
        richness = n_distinct(species),
        species = list(unique(species)),
        .groups = "drop"
      )
    
    if (!is.null(effort_by_day)) {
      out <- out %>% left_join(effort_by_day, by = "date")
    }
    
    if (return_species == TRUE) {
      species_table <- out %>%
        group_by(across(all_of(groups)), species) %>%
        summarise(.groups = "drop")
      
      return(list(richness = out, species = species_table))
    } else {
      return(out)
    }
  }
  
  # For seasonal richness (original behavior preserved, including join logic)
  if (per_temporal == "season") {
    if (is.null(spatial_group)) {
      out <- dt_filt %>%
        summarise(
          start_date = min(date),
          end_date   = max(date),
          richness   = n_distinct(species),
          species = list(unique(species))
        ) %>%
        left_join(effort_by_day, by = "date")
      
      if (return_species == TRUE) {
        species_table <- out %>%
          group_by(across(all_of(groups)), species) %>%
          summarise(.groups = "drop")
        
        return(list(richness = out, species = species_table))
      } else {
        return(out)
      }
    }
  }
}

run_spatial_rarefaction <- function(dt, B = 100, type = c("recorder", "site", "recorder_per_site")) {
  type <- match.arg(type)
  
  results <- list()
  
  if (type == "recorder") {
    recs <- unique(dt$recorder)
    results <- rbindlist(lapply(seq_along(recs), function(k) {
      rbindlist(lapply(seq_len(B), function(b) {
        sampled <- sample(recs, k)
        subset <- dt[recorder %in% sampled]
        r <- compute_richness(subset, schedule = NULL, per_spatial = "recorder", per_temporal = "day")
        alpha <- mean(r$richness)
        gamma <- length(unique(subset$species))
        data.table(
          unit = k,
          bootstrap = b,
          alpha = alpha,
          gamma = gamma,
          beta = gamma / alpha
        )
      }))
    }))
  } else if (type == "site") {
    sites <- unique(dt$site)
    results <- rbindlist(lapply(seq_along(sites), function(k) {
      rbindlist(lapply(seq_len(B), function(b) {
        sampled_sites <- sample(sites, k)
        subset <- dt[site %in% sampled_sites]
        r <- compute_richness(subset, schedule = NULL, per_spatial = "site", per_temporal = "day")
        alpha <- mean(r$richness)
        gamma <- length(unique(subset$species))
        data.table(
          unit = k,
          bootstrap = b,
          alpha = alpha,
          gamma = gamma,
          beta = gamma / alpha
        )
      }))
    }))
  } else if (type == "recorder_per_site") {
    
    results <- rbindlist(lapply(seq_len(B), function(b) {
      
      fractions <- seq(0.25, 1, by = 0.25)
      
      rbindlist(lapply(fractions, function(frac) {
        
        by_site <- split(dt, dt$site)
        richness_by_site <- rbindlist(lapply(by_site, function(site_dt) {
          
          recs <- unique(site_dt$recorder)
          n_to_sample <- max(1, floor(length(recs) * frac))
          
          sampled_recs <- sample(recs, n_to_sample)
          sub_dt <- site_dt[recorder %in% sampled_recs]
          
          r <- compute_richness(
            dt = sub_dt,
            schedule = NULL,
            per_spatial = "site",
            per_temporal = "day"
          )
          
          data.table(
            site = unique(site_dt$site),
            richness = mean(r$richness, na.rm = TRUE)
          )
        }))
        
        alpha <- mean(richness_by_site$richness, na.rm = TRUE)
        
        dt_all <- dt[recorder %in% unlist(lapply(by_site, function(s) {
          recs <- unique(s$recorder)
          n_to_sample <- max(1, floor(length(recs) * frac))
          sample(recs, n_to_sample)
        }))]
        
        gamma_dt <- compute_richness(
          dt = dt_all,
          schedule = NULL,
          per_spatial = "all",
          per_temporal = "day"
        )
        gamma <- mean(gamma_dt$richness, na.rm = TRUE)
        
        data.table(
          bootstrap = b,
          unit = round(frac * 100),
          alpha = alpha,
          gamma = gamma,
          beta = ifelse(alpha > 0, gamma / alpha, NA)
        )
      }))
    }))
    
    return(results)
  }
}

parse_path <- function(path, sep = "/") {
  parts <- str_split(path, sep)
  max_len <- max(lengths(parts))
  
  parts_mat <- t(sapply(parts, function(x) {
    length(x) <- max_len
    x
  }))
  
  as.data.frame(parts_mat, stringsAsFactors = FALSE)
}

split_path_matrix <- function(x, sep = "/") {
  if (length(x) == 0) return(NULL)
  parts <- stringr::str_split(x, sep)
  max_len <- max(lengths(parts))
  mat <- t(vapply(parts, function(p) {
    length(p) <- max_len
    p
  }, FUN.VALUE = character(max_len)))
  colnames(mat) <- paste0("V", seq_len(ncol(mat)))
  mat
}
split_path_matrix_preview <- function(x, sep = "/", n_preview = 20L) {
  if (length(x) == 0) return(NULL)
  x <- x[seq_len(min(length(x), n_preview))]
  split_path_matrix(x, sep = sep)
}
# ---- UI Definition ----
ui <- fluidPage(
  theme = shinytheme("darkly"),
  titlePanel("Acoustic sampling rarefaction"),
  tabsetPanel(
    # ---- TAB 1: data import----
    tabPanel(
      "1. Data import",
      sidebarLayout(
        sidebarPanel(
          h4("Import BirdNET data"),
          
          fileInput(
            inputId = "bn_file",
            label = "Upload BirdNET file(s) (.csv, .txt)",
            accept = c(".csv", ".txt"),
            multiple = TRUE
          ),
 
          
          # radioButtons("import_mode",
          #              "Import mode",
          #              c("File(s)" = "file",
          #                "Folder (recursive)" = "folder")),
          # 
          # conditionalPanel(
          #   "input.import_mode == 'file'",
          #   shinyFilesButton("bn_file", "Select BirdNET file(s)", "Choose", multiple = TRUE)
          # ),
          # 
          # conditionalPanel(
          #   "input.import_mode == 'folder'",
          #   shinyDirButton("bn_folder", "Select BirdNET folder", "Choose")
          # ),
          
          actionButton("load_data", "LOAD DATA AND PREVIW"),
          
          hr(),
          
          h4("Path parsing (optional)"),
          
          div(
            style = "overflow-x: auto; max-width: 100%;",
            tableOutput("beginpath_preview")
          ),
          
          uiOutput("path_parser_ui"),
          
          br(),
          
          actionButton("apply_parsing", "Parse path"),
          
          hr(),
          
          h4("Filter data"),
          
          numericInput("min_conf", "Min confidence", value = 0, min = 0, max = 1, step = 0.01),
          
          dateRangeInput(
            "date_subset",
            "Select time window",
            start = "2024-06-01",
            end = "2024-06-09"
          ),
          
          uiOutput("site_ui"),
          uiOutput("rec_ui"),
          uiOutput("year_ui"),
          uiOutput("species_ui"),
          
          actionButton("apply_filters", "Apply filters"),
          
          hr(),
          
          h4("Recorder coordinates"),
          
          radioButtons(
            "coord_mode",
            "Coordinate mode",
            choices = c(
              "Single coordinates (WGS84)" = "single",
              "Add multiple coordinates (file upload)" = "multiple"
            ),
            selected = "single"
          ),
          selectizeInput(
            "timezone",
            "Timezone ",
            choices = OlsonNames(),
            selected = "CET",
            options = list(
              placeholder = 'Type to search timezone...',
              maxOptions = 1000
            )
          ),
          
          conditionalPanel(
            condition = "input.coord_mode == 'single'",
            numericInput("single_lat",
                         "Latitude (WGS84)",
                         value = 46.19,
                         step = 0.0001),
            
            numericInput("single_lon",
                         "Longitude (WGS84)",
                         value = 8.13,
                         step = 0.0001)
          ),
          
          conditionalPanel(
            condition = "input.coord_mode == 'multiple'",
            fileInput("coord_file",
                      "Upload coordinate table (csv, txt, xlsx)",
                      accept = c(".csv", ".txt", ".xlsx")),
            helpText("File must contain columns: recorder OR site + latitude + longitude (WGS84)")
          )
        ),
        mainPanel(
          h4("Preview"),
          tableOutput("table_preview")
        )
      )
    ),
    
    # ---- TAB 2: Temporal Rarefaction ----
    tabPanel(
      "2. Duty cycle",
      fluidRow(
        column(
          width = 3,
          sidebarPanel(
            checkboxInput("record_24h", "Continuous 24h recording", FALSE),
            selectInput("start_event", "Start event", c("dawn", "sunrise")),
            selectInput("end_event", "End event", c("sunset", "dusk")),
            numericInput("length_morning", "Fixed duration (h)", 0, min = 0),
            numericInput("extend_before", "Extend before (h)", 1),
            numericInput("extend_after", "Extend after (h)", 0),
            checkboxInput("use_window2", "Add second daily window", FALSE), 
            conditionalPanel(
              condition = "input.use_window2 == true",
              selectInput("start_event2", "Start event (Window 2)", c("dawn", "sunrise", "sunset", "dusk")),
              selectInput("end_event2", "End event (Window 2)", c("dawn", "sunrise", "sunset", "dusk")),
              numericInput("length_morning2", "Fixed duration (h) - Window 2", 0, min = 0),
              numericInput("extend_before2", "Extend before (h) - Window 2", 0),
              numericInput("extend_after2", "Extend after (h) - Window 2", 0)), 
            numericInput("period_min", "Min period (min)", 1),
            numericInput("period_max", "Max period (min)", 30),
            numericInput("period_step", "Step (min)", 5),
            numericInput("duty_duration", "Duty duration (min)", 1),
            numericInput("B_temporal", "Bootstrap", value = 5, min = 1, step = 1),
            checkboxInput("show_summary", "Show average + CI", TRUE),
            actionButton("run_temporal", "Run analysis"),
            width = 12
          )
        ),
        column(
          width = 9,
          fluidRow(
            column(12, h4("Schedule preview")),
            column(6, plotOutput("sun_plot", height = "150px")),
            column(6, plotOutput("duty_plot", height = "150px"))
          ),
          br(),
          plotlyOutput("plot_temporal", height = "600px"),
          plotlyOutput("heatmap_missed_species", height = "600px")
        )
      )
    ),
    
    # ---- TAB 3: Time Window Rarefaction ----
    tabPanel(
      "3. Time window",
      fluidRow(
        column(
          width = 3,
          sidebarPanel(
            selectInput(
              "grid_event",
              "Reference event",
              choices = list(
                "Solar events" = c("dawn", "sunrise", "sunset", "dusk"),
                "Fixed hours"  = sprintf("%02d:00", 1:23)
              ),
              selected = "sunrise"
            ),
            
            numericInput("grid_start_min", "Min offset from event (hours)", value = -2, step = 0.25),
            numericInput("grid_start_max", "Max offset from event (hours)", value = 4, step = 0.25),
            numericInput("grid_start_step", "Start time step (hours)", value = 0.5, step = 0.25),
            numericInput("grid_dur_min", "Min duration (hours)", value = 0.5, step = 0.25),
            numericInput("grid_dur_max", "Max duration (hours)", value = 6, step = 0.25),
            numericInput("grid_dur_step", "Duration step (hours)", value = 0.5, step = 0.25),
            checkboxInput("use_window2_grid", "Add second daily window", FALSE), 
            
            conditionalPanel(
              condition = "input.use_window2_grid == true",
              
              selectInput(
                "grid_event2",
                "Reference event (Window 2)",
                choices = list(
                  "Solar events" = c("dawn", "sunrise", "sunset", "dusk"),
                  "Fixed hours"  = sprintf("%02d:00", 1:23)
                ),
                selected = "sunset"),
              numericInput("grid_start_min2", "Min offset (Window 2)", value = -2, step = 0.25),
              numericInput("grid_start_max2", "Max offset (Window 2)", value = 4, step = 0.25),
              numericInput("grid_start_step2", "Start time step (Window 2)", value = 0.5, step = 0.25),
              
              numericInput("grid_dur_min2", "Min duration (Window 2)", value = 0.5, step = 0.25),
              numericInput("grid_dur_max2", "Max duration (Window 2)", value = 6, step = 0.25),
              numericInput("grid_dur_step2", "Duration step (Window 2)", value = 0.5, step = 0.25)
              
            ), 
            actionButton("run_grid", "Run analysis"),
            width = 12,
            
            checkboxInput("show_frontier", "Show optimal frontier", value = TRUE)
          )
        ),
        column(
          width = 9,
          h4("Design preview"),
          uiOutput("grid_summary"),
          DT::DTOutput("grid_preview"),
          br(),
          h4("Richness heatmap"),
          uiOutput("grid_plots_ui")
          # plotOutput("grid_plot", height = "650px")
        )
      )
    ),
    
    # ---- TAB 4: multi optimum----
    tabPanel(
      "4. Multi optimum",
      fluidRow(
        column(
          width = 12,
          helpText(
            "This analysis combines the optimized duty cycle and daily window.",
            "Please run the Temporal and Window analyses first."
          ),
          actionButton("run_multi", "Run analysis")
        )
      ),
      fluidRow(
        column(
          width = 12,
          plotlyOutput("plot_multi_optimum", height = "650px")
        )
      )
    ),
    
    # ---- TAB 5: target species activity ----
    tabPanel(
      "5. Target species activity",
      sidebarLayout(
        sidebarPanel(
          width = 2,
          selectizeInput(
            "target_species",
            "Target species",
            choices = NULL,
            multiple = TRUE
          ),
          numericInput("target_period_min", "Period min", 1, min = 1),
          numericInput("target_period_max", "Period max", 60, min = 1),
          numericInput("target_period_step", "Period step", 5, min = 1),
          numericInput("target_duty_duration",
                       "Duty duration (min)",
                       value = 1,
                       min = 1),
          numericInput("target_nb_duty",
                       "Number of duty per period",
                       value = 1,
                       min = 1),
          selectInput(
            "time_resolution",
            "Time resolution plot",
            choices = c(5, 10, 15, 30, 45, 60),
            selected = 30
          ),
          actionButton("run_target", "Run analysis")
        ),
        mainPanel(
          width = 10,
          plotlyOutput("target_heatmap", height = "900px")
        )
      )
    ),
    
    # ---- TAB guidelines & help----
    tabPanel(
      "User Guide",
      
      fluidPage(
        
        h2("Application Overview"),
        p("whatever normal."),
        p(strong("BOLD.")),
        
        br(),
        
        h3("1. Importing Data"),
        p("coordinates multiple single, explantaion parsing path."),
        p(strong("Important:"), " parse =multiple coordinate, not parse = single coordinate."),
        
        br(),
        
        h3("2. Running Rarefaction"),
        p("explnanation bootstrap is a fixed nb of days (all season with replacmenet"),
        p(strong("Tip:"), " Start with small bootstrap for first look ."),
        
        br(),
        
        h3("3. window"),
        p("if you use window out or recording rabge, bias in result."),
        p(strong("Note:"), " whatever")
        
      )
    )
  )
)

# ---- Server Logic ----
server <- function(input, output, session) {
  message(Sys.time(), "session laucnched")
  # ---- values storage----
  values <- reactiveValues(
    summary_richness = NULL,
    grid_results     = NULL,
    schedules        = NULL,
    missed_summary   = NULL,
    experimental_domain = NULL
  )
  message(Sys.time(), "defining static data")
  # ---- recoprder battery information----
  battery_specs <- data.frame(
    recorder = c(
      "SMMicro 2",
      "SMMini 2",
      "SMMini 2",
      "SM4",
      "SM5",
      "AudioMoth -v1.2.0"
    ),
    nb_battery = c(4, 8, 6, 4, 11, 3),
    battery_type = c(
      "AA",
      "AA",
      "lithium",
      "4D alkaline",
      "lithium",
      "AA"
    ),
    duration_h = c(280, 530, 1330, 310, 1710, 189)
  )
  
  battery_specs$code <- paste0(
    gsub(" ", "", battery_specs$recorder),
    "_",
    battery_specs$nb_battery,
    substr(battery_specs$battery_type, 1, 2)
  )
  
  # ---- File selection roots ----
  message(Sys.time(), "root definition")
  # roots <- c(
  #   Home = normalizePath("~"),
  #   Documents = "~/Documents",
  #   Desktop = "~/Desktop"
  # )
  # message(Sys.time(), "choosing file")
  # shinyFileChoose(input, "bn_file", roots = roots)
  # shinyDirChoose(input, "bn_folder", roots = roots)
  
  # ---- Debug observer ----
  observeEvent(input$bn_folder, {
    cat("bn_folder changed:\n")
    str(input$bn_folder)
  })
  
  # ---- TAB 1: IMPORT/FORMAT ----
  bn_tables <- eventReactive(input$load_data, {
    
    req(input$bn_file)
    
    showNotification("Running: loading data...", type = "message", duration = 2)
    
    withProgress(message = "Loading and preparing BirdNET data...", value = 0, {
      
      message(Sys.time(), "starting load data")
      
      files <- input$bn_file$datapath
      
      incProgress(0.2, detail = "Reading files...")
      
      # data_list <- lapply(files, read.delim)
      data_list <- lapply(files, fread) 
      incProgress(0.6, detail = "Combining files...")
      result <- rbindlist(data_list, fill = TRUE)
      message(Sys.time(), "standradize column names")
        names(result) <- make.names(names(result))
        message(Sys.time(), "standradize column names - done")
        # Original normalization
        if ("Begin.Path" %in% names(result)) {
          message(Sys.time(), "standardize begin path")
          result[, Begin.Path := normalize_begin_path(Begin.Path)]
        }
        # incProgress(0.2, detail = "Formatting dataset 1/2...")
        result

      incProgress(1)
      
      result
      
    })
  })
  # bn_tables <- eventReactive(input$load_data, {
  #   showNotification("Running: loading data...", type="message", duration=2)
  #   withProgress(message = "Loading and preparing BirdNET data...", value = 0, {
  #     message(Sys.time(), "starting load data")
  #     files <- NULL
  #     
  #     if (input$import_mode == "file") {
  #       message(Sys.time(), "file mode")
  #       paths <- parseFilePaths(roots, input$bn_file)
  #       if (!is.null(paths) && nrow(paths) > 0) files <- paths$datapath
  #       
  #     } else if (input$import_mode == "folder") {
  #       message(Sys.time(), "folder mode")
  #       folder <- parseDirPath(roots, input$bn_folder)
  #       
  #       files <- list.files(
  #         folder,
  #         recursive = TRUE,
  #         full.names = TRUE,
  #         pattern = "\\.selection\\.table\\.txt$"
  #       )
  #     }
  #     
      # validate(need(length(files) > 0, "No BirdNET files found"))
      
      # withProgress(message = "Reading BirdNET files...", value = 0, {
    #   message(Sys.time(), "binding multiples file")
    #   n_files <- length(files)
    #   tables <- lapply(seq_along(files), function(i) {
    #     incProgress(1 / n_files, detail = paste(i, "/", n_files))
    #     fread(files[i])
    #   })
    #   
    #   incProgress(0.3, detail = "Binding tables...")
    #   result <- rbindlist(tables, fill = TRUE)
    #   
    #   incProgress(0.5, detail = "Standardizing column names...")
    #   message(Sys.time(), "standradize column names")
    #   names(result) <- make.names(names(result))
    #   message(Sys.time(), "standradize column names - done")
    #   # Original normalization
    #   if ("Begin.Path" %in% names(result)) {
    #     message(Sys.time(), "standardize begin path")
    #     result[, Begin.Path := normalize_begin_path(Begin.Path)]
    #   }
    #   incProgress(0.2, detail = "Formatting dataset 1/2...")
    #   result
    # })
# })  
  
  observeEvent(input$load_data, {
    updateSelectInput(session, "site_col", selected = "None")
    updateSelectInput(session, "rec_col", selected = "None")
  })
  
  # ---- extracting metadata ----
  dt_parsed <- reactive({
    req(bn_tables())
    message(Sys.time(), Sys.time(), "starting to format data and columns")
    withProgress(message = "Loading and preparing BirdNET data...", value = 0, {
      # withProgress(message = "Formatting BirdNET data...", value = 0, {
      
      dt <- bn_tables()
      
      # Defensive checks (SAFE: errors earlier with clear message)
      need_cols(dt, c("Begin.Path", "Common.Name", "Confidence", "Begin.Time..s."), context = "BirdNET selection table")
      message(Sys.time(), Sys.time(), "remove extension")
      # Extract clean filename (no extension)
      dt[, filename := tools::file_path_sans_ext(basename(Begin.Path))]
      
      incProgress(0.2, detail = "Extracting metadata")
      message(Sys.time(), Sys.time(), "big formatting starting")
      # Original transformation logic preserved
      BIRDNET <- dt %>%
        mutate(
          species = Common.Name,
          conf = Confidence,
          start = Begin.Time..s.,
          timestamp = str_extract(filename, "\\d{8}_\\d{6}"),
          date = paste0(substr(timestamp, 1, 4), "-", substr(timestamp, 5, 6), "-", substr(timestamp, 7, 8)),
          time = paste0(substr(timestamp, 10, 11), ":", substr(timestamp, 12, 13), ":", substr(timestamp, 14, 15)),
          # date = str_replace(timestamp,
          #                    "^(\\d{4})(\\d{2})(\\d{2})_(\\d{2})(\\d{2})(\\d{2})$",
          #                    "\\1-\\2-\\3"),
          # time = str_replace(timestamp,
          #                    "^(\\d{4})(\\d{2})(\\d{2})_(\\d{2})(\\d{2})(\\d{2})$",
          #                    "\\4:\\5:\\6"),
          timestamp_fm = ymd_hms(timestamp, tz = "UTC"),
          timestamp_adjusted = timestamp_fm + seconds(start),
          timestamp = format(timestamp_adjusted, "%Y-%m-%d %H:%M:%S")
        ) %>%
        filter(!is.na(timestamp)) %>%
        mutate(
          year = year(date),
          minute_of_day = hour(timestamp_adjusted) * 60 + minute(timestamp_adjusted),
          start_file = time,
          site = "all",
          recorder = "global"
        ) %>%
        select(species, site, recorder, filename, start,
               conf, date, year, start_file, timestamp, minute_of_day, Begin.Path)
      
      
      BIRDNET <- as.data.table(BIRDNET)
      message(Sys.time(), "big formatting - done ")
      incProgress(0.8, detail = "finalizing")
      BIRDNET
    })
  })
  
  # ---- Begin.Path preview and UI----
  output$beginpath_preview <- renderTable({
    
    req(bn_tables())
    message(Sys.time(), "begin path preview")
    dt <- bn_tables()
    if (!"Begin.Path" %in% names(dt)) return(NULL)
    message(Sys.time(), "splitting path")
    incProgress(0.8, detail = "splitting path preview")
    mat <- split_path_matrix_preview(dt$Begin.Path)
    preview <- as.data.frame(mat)
    message(Sys.time(), "building preview")
    colnames(preview) <- paste0("V", seq_len(ncol(preview)))
    
    head(preview, 6)
  })
  
  # ---- Path parser UI ----
  output$path_parser_ui <- renderUI({
    
    req(bn_tables())
    
    dt <- bn_tables()
    message(Sys.time(), "building path parser ")
    if (!"Begin.Path" %in% names(dt)) return(tags$div("Begin Path column missing"))
    
    mat <- split_path_matrix_preview(dt$Begin.Path)
    n <- ncol(mat)
    message(Sys.time(), "building path parser - choices")
    tagList(
      selectInput("site_col", "Site column",
                  c("None", paste0("V", seq_len(n)))),
      selectInput("rec_col", "Recorder column",
                  c("None", paste0("V", seq_len(n))))
    )
    # message(Sys.time(), "building path parser  - done ")
    # })
  })
  
  # ---- apply parsing ----
  dt_parsed_path <- reactiveVal(NULL)
  preview_data <- reactiveVal(NULL)
  observeEvent(input$apply_parsing, {
    
    req(dt_parsed())
    message(Sys.time(), "aplliying path parsing")
    withProgress(message = "Applying path parsing...", value = 0, {
      
      dt <- copy(dt_parsed())
      message(Sys.time(), "aplliying path parsing, split 2")
      mat <- split_path_matrix(dt$Begin.Path)
      
      site_col <- if (!is.null(input$site_col) && input$site_col != "None")
        as.integer(sub("V", "", input$site_col))
      else NULL
      
      rec_col <- if (!is.null(input$rec_col) && input$rec_col != "None")
        as.integer(sub("V", "", input$rec_col))
      else NULL
      
      incProgress(0.5)
      message(Sys.time(), "aplliying path parsing, np parse global all")
      dt[, site := if (!is.null(site_col)) mat[, site_col] else "all"]
      dt[, recorder := if (!is.null(rec_col)) mat[, rec_col] else "global"]
      
      incProgress(0.5)
      
      dt_parsed_path(dt)
      
      print(unique(dt$site))
      
      print(unique(dt$recorder))
    })
  })
  
  observeEvent(input$load_data, {
    dt_parsed_path(NULL)
  })
  
  # ---- reload dataset if updated----
  dt_for_app <- reactive({
    message(Sys.time(), "replacing parsed dataset if reloaded")
    req(dt_parsed())
    
    parsed <- dt_parsed_path()
    if (!is.null(parsed)) return(parsed)
    
    dt_parsed()
  })
  dt_filtered <- reactiveVal(NULL)
  
  # ---- Reset filters dataset ----
  observeEvent(input$load_data, {
    message(Sys.time(), "reset filtered dataset on new load")
    dt_filtered(NULL)
  })
  
  dt_raw <- reactive({
    req(dt_for_app())
    
    # If no filters applied yet -> return full dataset for preview/analysis
    filtered <- dt_filtered()
    if (is.null(filtered)) {
      return(dt_for_app())
    }
    
    filtered
  })
  # ---- filtering data ----
  
  observeEvent(input$apply_filters, {
    req(dt_for_app())
    message(Sys.time(), "apllying filters")
    
    withProgress(message = "Filtering data...", value = 0, {
      
      dt <- copy(dt_for_app())
      
      incProgress(0.2, detail = "Filtering by date")
      message(Sys.time(), "apllying filters- date")
      if (!is.null(input$date_subset) && length(input$date_subset) == 2) {
        dt <- dt[
          as.Date(date) >= input$date_subset[1] &
            as.Date(date) <= input$date_subset[2]
        ]
      }
      incProgress(0.1, detail = "Filtering by confidence")
      if (!is.null(input$min_conf)) dt <- dt[conf >= input$min_conf]
      message(Sys.time(), "apllying filters- site/recorder")
      incProgress(0.1, detail = "Filtering by site / recorder")
      if (!is.null(input$site)) dt <- dt[site %in% input$site]
      if (!is.null(input$recorder)) dt <- dt[recorder %in% input$recorder]
      
      incProgress(0.2, detail = "Filtering by year")
      message(Sys.time(), "apllying filters- year")
      if (!is.null(input$year)) dt <- dt[year %in% input$year]
      
      message(Sys.time(), "apllying filters- species")
      incProgress(0.2, detail = "Filtering by species")
      if (!is.null(input$species)) dt <- dt[species %in% input$species]
      
      message(Sys.time(), "apllying filters- done ")
      incProgress(0.2, detail = "Done")
      
      # store filtered dataset (this is the whole point)
      dt_filtered(dt)
    })
  })
  
  
  # ---- Filter selectors ----
  output$site_ui <- renderUI({
    req(dt_for_app())
    sites <- unique(dt_for_app()$site)
    selectInput("site", "Site", sites, multiple = TRUE)
  })
  
  output$rec_ui <- renderUI({
    req(dt_for_app())
    selectInput("recorder", "Recorder", unique(dt_for_app()$recorder), multiple = TRUE)
  })
  
  output$year_ui <- renderUI({
    req(dt_for_app())
    selectInput("year", "Year", unique(dt_for_app()$year), multiple = TRUE)
  })
  
  output$species_ui <- renderUI({
    req(dt_for_app())
    selectInput("species", "Species", sort(unique(dt_for_app()$species)), multiple = TRUE)
  })
  
  observeEvent(input$load_data, {
    message(Sys.time(), "update selection site and rec")
    if (!is.null(input$site_col)) updateSelectInput(session, "site_col", selected = "None")
    if (!is.null(input$rec_col))  updateSelectInput(session, "rec_col",  selected = "None")
    message(Sys.time(), "aplly all global to rec and site if null")
    if (!is.null(input$site))     updateSelectInput(session, "site",     selected = "all")
    if (!is.null(input$recorder)) updateSelectInput(session, "recorder", selected = "global")
  })
  
  # ---- display final preview ----
  output$table_preview <- renderTable({
    
    showNotification("time is long 7", type = "message", duration = 5)
    
    req(dt_raw())
    head(dt_raw(), 20)
  })
  
  
  # ---- define coordinates----
  coord_single <- reactive({
    req(input$coord_mode == "single")
    message(Sys.time(), "coord mode single")
    data.frame(latitude = input$single_lat, longitude = input$single_lon)
  })
  
  coord_multiple <- reactive({
    req(input$coord_mode == "multiple")
    req(input$coord_file)
    message(Sys.time(), "coord mode multiple")
    df <- read_coord_file(input$coord_file$datapath, input$coord_file$name)
    df
  })
  
  # ---- Active coordinates ----
  active_lat <- reactive({
    if (input$coord_mode == "single") {
      input$single_lat
    } else {
      req(coord_multiple())
      coord_multiple()$latitude[1]
    }
  })
  
  active_lon <- reactive({
    if (input$coord_mode == "single") {
      input$single_lon
    } else {
      req(coord_multiple())
      coord_multiple()$longitude[1]
    }
  })
  
  # ---- pre-compute solar table----
  solar_table <- reactive({
    
    req(dt_for_app())
    message(Sys.time(), "building solar table from coord")
    req(input$timezone)
    
    dt <- dt_for_app()
    dates_unique <- unique(dt$date)
    
    # CASE 1 : SINGLE COORDINATE
    if (input$coord_mode == "single") {
      message(Sys.time(), "sun table for single")
      lat <- input$single_lat
      lon <- input$single_lon
      
      solar_list <- lapply(dates_unique, function(d) {
        
        sun <- getSunlightTimes(
          date = as.Date(d),
          lat = lat,
          lon = lon,
          keep = c("dawn", "sunrise", "sunset", "dusk"),
          tz = input$timezone
        )
        
        data.table(
          date = as.Date(d),
          recorder = "global",
          dawn = sun$dawn,
          sunrise = sun$sunrise,
          sunset = sun$sunset,
          dusk = sun$dusk
        )
      })
      message(Sys.time(), "sun table for single - done")
      return(data.table::rbindlist(solar_list))
    }
    
    # CASE 2 : MULTIPLE COORDINATES
    if (input$coord_mode == "multiple") {
      
      req(coord_multiple())
      
      coords <- coord_multiple()
      
      req("recorder" %in% names(coords))
      req("latitude" %in% names(coords))
      req("longitude" %in% names(coords))
      message(Sys.time(), "sun table for multiple")
      comb <- unique(dt[, .(recorder, date)])
      comb <- merge(comb, coords, by = "recorder", all.x = TRUE)
      
      solar_list <- lapply(seq_len(nrow(comb)), function(i) {
        
        sun <- getSunlightTimes(
          date = as.Date(comb$date[i]),
          lat = comb$latitude[i],
          lon = comb$longitude[i],
          keep = c("dawn", "sunrise", "sunset", "dusk"),
          tz = input$timezone
        )
        
        data.table(
          recorder = comb$recorder[i],
          date = as.Date(comb$date[i]),
          dawn = sun$dawn,
          sunrise = sun$sunrise,
          sunset = sun$sunset,
          dusk = sun$dusk
        )
      })
      message(Sys.time(), "sun table for multiple - done")
      return(data.table::rbindlist(solar_list))
    }
  })
  
  # ---- TAB 2: Temporal rarefaction ----
  # ---- sun & duty plot ----
  output$sun_plot <- renderPlot({
    req(input$start_event, input$end_event, input$extend_before, input$extend_after)
    message(Sys.time(), "sun schedule previw")
    showNotification("Schedule preview updated", type = "message", duration = 1)
    
    sun_events <- getSunlightTimes(
      date = Sys.Date(),
      lat = active_lat(),
      lon = active_lon(),
      keep = c("dawn", "sunrise", "sunset", "dusk"),
      tz = input$timezone
    )
    
    sun_events <- as.data.table(sun_events)
    sun_events <- melt(sun_events, measure.vars = c("dawn", "sunrise", "sunset", "dusk"),
                       variable.name = "event", value.name = "datetime")
    setDT(sun_events)
    sun_events[, minute_of_day := hour(datetime) * 60 + minute(datetime)]
    
    
    
    # ---- Window 1 ----
    start_min1 <- sun_events[event == input$start_event, minute_of_day]
    
    if (input$length_morning > 0) {
      end_min1 <- start_min1 + input$length_morning * 60
    } else {
      end_min1 <- sun_events[event == input$end_event, minute_of_day]
    }
    
    start_min1 <- start_min1 - input$extend_before * 60
    end_min1   <- end_min1   + input$extend_after  * 60
    #changement de jour
    start_min1 <- start_min1 %% (24*60)
    end_min1   <- end_min1   %% (24*60)
    
    # ---- 24H mode ----
    if (isTRUE(input$record_24h)) {
      
      start_min1 <- 0
      end_min1   <- 1440
      
    }
    # ---- 2 window case ----
    if (isTRUE(input$use_window2)) {
      
      start_min2 <- sun_events[event == input$start_event2, minute_of_day]
      
      if (input$length_morning2 > 0) {
        end_min2 <- start_min2 + input$length_morning2 * 60
      } else {
        end_min2 <- sun_events[event == input$end_event2, minute_of_day]
      }
      
      start_min2 <- start_min2 - input$extend_before2 * 60
      end_min2   <- end_min2   + input$extend_after2  * 60
      start_min2 <- start_min2 %% (24*60)
      end_min2   <- end_min2   %% (24*60)
    }
    
    message(Sys.time(), "sun schedule preview calculation done ")
    
    
    seg1 <- split_window_across_days(Sys.Date(), start_min1, end_min1)
    print(start_min1)
    print(end_min1)
    p <- ggplot()
    
    # ---- Sun plot----
    if (start_min1 <= end_min1) {
      
      p <- p +
        geom_rect(
          aes(xmin = start_min1, xmax = end_min1),
          ymin = 0.9, ymax = 1.1,
          fill = "skyblue", alpha = 0.5
        )
      
    } else {
      
      p <- p +
        geom_rect(
          aes(xmin = start_min1, xmax = 1440),
          ymin = 0.9, ymax = 1.1,
          fill = "skyblue", alpha = 0.5
        ) +
        geom_rect(
          aes(xmin = 0, xmax = end_min1),
          ymin = 0.9, ymax = 1.1,
          fill = "skyblue", alpha = 0.5
        )
    }
    
    # ---- Sun plot (Window 2 )----
    if (isTRUE(input$use_window2)) {
      
      if (start_min2 <= end_min2) {
        
        p <- p +
          geom_rect(
            aes(xmin = start_min2, xmax = end_min2),
            ymin = 0.9, ymax = 1.1,
            fill = "orange", alpha = 0.5
          )
        
      } else {
        
        p <- p +
          geom_rect(
            aes(xmin = start_min2, xmax = 1440),
            ymin = 0.9, ymax = 1.1,
            fill = "orange", alpha = 0.5
          ) +
          geom_rect(
            aes(xmin = 0, xmax = end_min2),
            ymin = 0.9, ymax = 1.1,
            fill = "orange", alpha = 0.5
          )
      }
    }
    p +
      geom_vline(data = sun_events,
                 aes(xintercept = minute_of_day, color = event),
                 linewidth = 0.8) +
      scale_color_manual(values = c(
        "dawn" = "darkorange",
        "sunrise" = "gold",
        "sunset" = "firebrick",
        "dusk" = "purple"
      )) +
      scale_x_continuous(limits = c(0, 1440), breaks = NULL)+
      # scale_x_continuous(limits = c(200, 1300), breaks = NULL) +
      theme_minimal() +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank()
      )
  })
  
  
  # ---- Duty plot----
  output$duty_plot <- renderPlot({
    req(input$period_min, input$period_max, input$period_step, input$duty_duration)
    message(Sys.time(), "sun duty cycle preview")
    showNotification("Duty cycle plot updated", type = "message", duration = 1)
    
    # start_min <- 360
    # end_min   <- 1080
    start_min <- 0
    end_min   <- 24*60
    
    periods <- seq(input$period_min, input$period_max, by = input$period_step)
    duty_duration <- input$duty_duration
    
    all_duty_cycles <- rbindlist(
      lapply(periods, function(p) {
        start_times <- seq(start_min, end_min - duty_duration, by = p)
        data.table(
          period = p,
          start = start_times,
          end = start_times + duty_duration
        )
      })
    )
    
    all_duty_cycles[, y := factor(period, levels = rev(unique(period)))]
    message(Sys.time(), "sun duty cycle preview, calculation done ")
    ggplot(all_duty_cycles) +
      geom_rect(
        aes(
          xmin = start, xmax = end,
          ymin = as.numeric(y) - 0.4,
          ymax = as.numeric(y) + 0.4
        ),
        fill = "skyblue", alpha = 1
      ) +
      scale_y_continuous(
        breaks = seq_along(rev(periods)),
        labels = paste0("Period = ", rev(periods), " min")
      ) +
      theme_void() +
      theme(
        axis.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()
      )
  })
  
  # ---- TAB 2 :COMPUTATION ----
  
  temporal_results <- eventReactive(input$run_temporal, {
    
    req(input$period_step > 0)
    message(Sys.time(), "temporal raref start")
    # ---- creating schedules----
    periods <- seq(
      as.numeric(input$period_min),
      as.numeric(input$period_max),
      by = as.numeric(input$period_step)
    )
    message(Sys.time(), "period sequence started")
    param_grid <- CJ(period = periods, duty_duration = input$duty_duration)
    message(Sys.time(), "grid built")
    local_dt <- dt_raw()
    req(local_dt)
    message(Sys.time(), "temporal raref start - creating schedules")
    withProgress(message = "Creating schedules...", value = 0, {
      
      recs <- unique(local_dt$recorder)
      solar_dt <- solar_table()
      solar_dt <- setDT(solar_dt)
      
      schedules <- lapply(seq_len(nrow(param_grid)), function(i) {
        
        p <- param_grid[i]
        
        sched_rec_list <- lapply(recs, function(r) {
          
          sun_r <- solar_dt[solar_dt$recorder == r, ]
          
          sched_r1 <- create_sun_duty_schedule(
            sun_table = sun_r,
            start_event = input$start_event,
            end_event   = input$end_event,
            length_morning = input$length_morning,
            extend_before_hours = input$extend_before,
            extend_after_hours  = input$extend_after,
            period = p$period,
            duty_duration = p$duty_duration,
            nb_duty = 1,
            random_start = FALSE,
            record_24h = input$record_24h
          )
          if (isTRUE(input$use_window2)) {
            
            sched_r2 <- create_sun_duty_schedule(
              sun_table = sun_r,
              start_event = input$start_event2,
              end_event   = input$end_event2,
              length_morning = input$length_morning2,
              extend_before_hours = input$extend_before2,
              extend_after_hours  = input$extend_after2,
              period = p$period,
              duty_duration = p$duty_duration,
              nb_duty = 1,
              random_start = FALSE
            )
            
            sched_r <- rbindlist(list(sched_r1, sched_r2))
            
          } else {
            
            sched_r <- sched_r1
          }
          sched_r[, recorder := r]
          sched_r
        })
        
        sched <- rbindlist(sched_rec_list)
        setattr(sched, "period", p$period)
        setattr(sched, "duty_duration", p$duty_duration)
        
        sched
      })
      message(Sys.time(), "schedules stored")
      values$schedules <- schedules
      
      
      schedules
    })
    message(Sys.time(), "schedules done")
    # ---- store schedule----
    values$experimental_domain <- list(
      
      timezone = input$timezone,
      
      w1 = list(
        start_event   = input$start_event,
        end_event     = input$end_event,
        length        = input$length_morning,
        extend_before = input$extend_before,
        extend_after  = input$extend_after
      ),
      
      w2 = if (isTRUE(input$use_window2)) list(
        start_event   = input$start_event2,
        end_event     = input$end_event2,
        length        = input$length_morning2,
        extend_before = input$extend_before2,
        extend_after  = input$extend_after2
      ) else NULL
    )
    
    
    # ---- bootstrap comoutation----
    message(Sys.time(), "experimental design stored ")
    withProgress(message = "Running bootstrap iterations...", value = 0, {
      
      all_results <- vector("list", input$B_temporal)
      message(Sys.time(), "temporal raref - bootstraping")
      for (b in seq_len(input$B_temporal)) {
        
        incProgress(1 / input$B_temporal, detail = paste(b, "/", input$B_temporal))
        
        dates <- unique(local_dt$date)
        sampled_days <- sample(dates, length(dates), replace = TRUE)
        dt_boot <- local_dt[local_dt$date %in% sampled_days]
        w0 <- values$schedules[[1L]][, min(start_min)]
        w1 <- values$schedules[[1L]][, max(end_min)]
        
        if (w0 <= w1) {
          
          # fenêtre normale
          dt_boot <- dt_boot[
            minute_of_day >= w0 &
              minute_of_day <= w1
          ]
          
        } else {
          
          # fenêtre traverse minuit
          dt_boot <- dt_boot[
            minute_of_day >= w0 |
              minute_of_day <= w1
          ]
        }
        # w0 <- values$schedules[[1L]][, min(start_min)]
        # w1 <- values$schedules[[1L]][, max(end_min)]
        # dt_boot <- dt_boot[minute_of_day >= w0 & minute_of_day <= w1]
        
        schedule_ref <- values$schedules[[1]]
        
        ref_result <- compute_richness(
          dt = dt_boot,
          schedule = schedule_ref,
          per_spatial = "all",
          per_temporal = "day",
          return_species = TRUE
        )
        
        ref_species <- unique(ref_result$species$species)
        
        one_bootstrap <- lapply(values$schedules, function(s) {
          
          res <- compute_richness(
            dt = dt_boot,
            schedule = s,
            per_spatial = "all",
            per_temporal = "day",
            return_species = TRUE
          )
          
          r <- res$richness
          sp <- res$species
          
          r_boot <- r[r$date %in% sampled_days, ]
          sp_boot <- sp[sp$date %in% sampled_days, ]
          
          data.table(
            period = attr(s, "period"),
            effort = s[, sum(end_min - start_min + 1)],
            richness = mean(r_boot$richness),
            bootstrap = b,
            species_detected = list(unique(sp_boot$species))
          )
        })
        
        all_results[[b]] <- rbindlist(one_bootstrap)
      }
      message(Sys.time(), "temp raref loop done ")
      summary_richness <- rbindlist(all_results)
      
      summary_richness[, richness_pct :=
                         100 * richness / max(richness),
                       by = bootstrap
      ]
      
      # ---- missed species list---
      message(Sys.time(), "temp raref - relative richensss claculated")
      reference_table <- summary_richness[
        period == min(period),
        .(reference_species = list(unlist(species_detected))),
        by = bootstrap
      ]
      message(Sys.time(), "temp raref missing species ")
      summary_richness <- merge(
        summary_richness,
        reference_table,
        by = "bootstrap",
        all.x = TRUE
      )
      message(Sys.time(), "temp raref missing species 2 ")
      missed_species_list <- summary_richness[, .(
        species_missed = list(setdiff(reference_species[[1]], unlist(species_detected)))
      ), by = .(bootstrap, period)]
      message(Sys.time(), "temp raref missing species 3")
      missed_species_long <- missed_species_list[, .(
        species_missed = unlist(species_missed)
      ), by = .(bootstrap, period)]
      message(Sys.time(), "temp raref missing species 4")
      missed_summary <- missed_species_long[, .N, by = .(species_missed, period)]
      message(Sys.time(), "temp raref missing species 5 ")
      species_order <- missed_summary[
        , .(total_missed = sum(N)), by = species_missed
      ][order(total_missed), species_missed]
      message(Sys.time(), "temp raref missing species 6")
      missed_summary[, species_missed := factor(species_missed, levels = species_order)]
      message(Sys.time(), "temp raref missing species 7 ")
      values$missed_summary <- missed_summary
      values$summary_richness <- summary_richness
      message(Sys.time(), "temp raref calculation - done  ")
      summary_richness
    })
  })
  
  # ---- TAB 2 : PLOT TEMPORAL RAREF ----
  output$plot_temporal <- renderPlotly({
    req(temporal_results())
    
    validate(
      need(
        nrow(temporal_results()) > 0,
        "Invalid temporal window.\nEnd time occurs before start time.\nCheck start/end events or extensions."
      )
    )
    message(Sys.time(), "temp raref plot")
    ci_summary <- temporal_results()[, .(
      mean_richness = mean(richness_pct, na.rm = TRUE),
      lower = quantile(richness_pct, probs = 0.025, na.rm = TRUE),
      upper = quantile(richness_pct, probs = 0.95, na.rm = TRUE)
    ), by = effort]
    
    durations <- battery_specs$duration_h
    codes     <- battery_specs$code
    # ---- plot rarefaction----
    gg <- ggplot(
      temporal_results(),
      aes(x = effort, y = richness_pct, group = bootstrap, color = factor(bootstrap))
    ) +
      geom_line(alpha = 0.5) +
      geom_point(alpha = 0.5) +
      labs(x = "sampling schedule - effort (min)", y = "relative richness [%]") +
      theme_minimal() +
      guides(color = "none")
    
    if (input$show_summary) {
      gg <- gg +
        geom_ribbon(
          data = ci_summary,
          mapping = aes(x = .data$effort, ymin = .data$lower, ymax = .data$upper),
          inherit.aes = FALSE,
          fill = "skyblue", alpha = 0.3
        ) +
        geom_line(
          data = ci_summary,
          mapping = aes(
            x = .data$effort,
            y = .data$mean_richness,
            text = paste0(
              "Effort: ", .data$effort, " min - ", round(.data$effort / 60, 2), " hours ",
              "<br>Richness: ", round(.data$mean_richness, 1), "%",
              "<br><br><b>Full battery sets needed:</b><br>",
              paste0(
                codes[1], ": ", round((.data$effort / 60) / durations[1], 2), "<br>",
                codes[2], ": ", round((.data$effort / 60) / durations[2], 2), "<br>",
                codes[3], ": ", round((.data$effort / 60) / durations[3], 2), "<br>",
                codes[4], ": ", round((.data$effort / 60) / durations[4], 2), "<br>",
                codes[5], ": ", round((.data$effort / 60) / durations[5], 2), "<br>",
                codes[6], ": ", round((.data$effort / 60) / durations[6], 2)
              )
            )
          ),
          inherit.aes = FALSE,
          color = "blue", linewidth = 1.2
        )
    }
    
    ggplotly(gg, tooltip = "text")
  })
  
  message(Sys.time(), "temp raref curve done")
  # ---- plot missed species----
  output$heatmap_missed_species <- renderPlotly({
    
    req(values$missed_summary)
    message(Sys.time(), "temp raref missing species plot ")
    missed_summary <- values$missed_summary
    period_min <- min(missed_summary$period)
    period_max <- max(missed_summary$period)
    period_step <- unique(diff(sort(unique(missed_summary$period))))[1]
    message(Sys.time(), "temp raref missing species plot 2  ")
    mm <- ggplot(missed_summary, aes(x = period, y = species_missed, fill = N)) +
      geom_tile(aes(text = paste0(
        "Species: ", species_missed,
        "\nPeriod: ", period, " min",
        "\nTimes missed: ", N
      )), color = "white", width = period_step, height = 1) +
      scale_fill_viridis_c(option = "D", name = "# bootstraps missed", na.value = "white") +
      scale_x_continuous(
        name = "Period (min)",
        breaks = seq(period_min, period_max, by = period_step),
        expand = c(0, 0)
      ) +
      labs(
        y = "Species",
        title = "Species missed across duty cycles"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.major = element_line(color = "grey85", linewidth  = 0.3),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank()
      ) +
      geom_vline(
        xintercept = seq(
          period_min - period_step / 2,
          period_max + period_step / 2,
          by = period_step
        ),
        color = "grey90",
        linewidth = 0.2
      ) +
      geom_hline(yintercept = seq(0.5, length(levels(missed_summary$species_missed)) + 0.5, by = 1),
                 color = "grey90", linewidth = 0.2)
    
    ggplotly(mm, tooltip = "text") %>%
      plotly::config(
        modeBarButtonsToRemove = c("select2d", "lasso2d", "zoomIn2d", "zoomOut2d"),
        displaylogo = FALSE
      )
  })
  
  message(Sys.time(), "temp raref missing species  plot - done ")
  
  # ---- TAB 3: DAILY WINDOW----
  
  output$grid_preview <- DT::renderDataTable({
    # ---- grid preview--
    
    req(input$grid_start_min,
        input$grid_start_max,
        input$grid_start_step,
        input$grid_dur_min,
        input$grid_dur_max,
        input$grid_dur_step)
    
    
    # ---- Window 1 grid ----
    starts <- seq(input$grid_start_min,
                  input$grid_start_max,
                  by = input$grid_start_step)
    
    durs <- seq(input$grid_dur_min,
                input$grid_dur_max,
                by = input$grid_dur_step)
    
    design <- as.data.table(expand.grid(
      start_offset_h = starts,
      duration_h     = durs
    ))
    
    design[, end_offset_h := start_offset_h + duration_h]
    
    design[, Window := sprintf("%+.1fh → %+.1fh",
                               start_offset_h,
                               end_offset_h)]
    
    design[, Effort := sprintf("%.1f h", duration_h)]
    
    design[, Start := sprintf("Start = Reference event %+.1f h",
                              start_offset_h)]
    
    
    # ---- Window 2 grid (optional) ----
    if (isTRUE(input$use_window2_grid)) {
      
      req(input$grid_start_min2,
          input$grid_start_max2,
          input$grid_start_step2,
          input$grid_dur_min2,
          input$grid_dur_max2,
          input$grid_dur_step2)
      
      starts2 <- seq(input$grid_start_min2,
                     input$grid_start_max2,
                     by = input$grid_start_step2)
      
      durs2 <- seq(input$grid_dur_min2,
                   input$grid_dur_max2,
                   by = input$grid_dur_step2)
      
      design2 <- as.data.table(expand.grid(
        start_offset_h2 = starts2,
        duration_h2     = durs2
      ))
      
      design2[, end_offset_h2 := start_offset_h2 + duration_h2]
      
      design2[, Window2 := sprintf("%+.1fh → %+.1fh",
                                   start_offset_h2,
                                   end_offset_h2)]
      
      design2[, Effort2 := sprintf("%.1f h", duration_h2)]
      
      
      design <- cbind(design, design2[1:nrow(design)])
    }
    
    setorder(design, start_offset_h, duration_h)
    
    # ---- grid preview UI ----
    if (isTRUE(input$use_window2_grid)) {
      
      out <- design[, .(
        start_group = Start,
        Window,
        Effort,
        Window2,
        Effort2
      )]
      
    } else {
      
      out <- design[, .(
        start_group = Start,
        Window,
        Effort
      )]
    }
    
    DT::datatable(
      out,
      rownames = FALSE,
      extensions = "RowGroup",
      options = list(
        paging = FALSE,
        scrollY = "150px",
        scrollCollapse = TRUE,
        info = FALSE,
        searching = FALSE,
        ordering = FALSE,
        dom = "t",
        rowGroup = list(dataSrc = 0),
        columnDefs = list(
          list(visible = FALSE, targets = 0),
          list(className = "dt-center", targets = "_all")
        ),
        class = "compact cell-border"
      )
    ) %>%
      formatStyle(names(out), color = "white")
    
  })
  
  # ---- text summary for table----
  output$grid_summary <- renderUI({
    
    req(input$grid_start_min, input$grid_start_max, input$grid_start_step,
        input$grid_dur_min, input$grid_dur_max, input$grid_dur_step)
    message(Sys.time(), "window tab - grid UI")
    n_starts <- length(seq(input$grid_start_min, input$grid_start_max, by = input$grid_start_step))
    n_durs <- length(seq(input$grid_dur_min, input$grid_dur_max, by = input$grid_dur_step))
    if (isTRUE(input$use_window2_grid)) {
      
      n_starts2 <- length(seq(input$grid_start_min2,
                              input$grid_start_max2,
                              by = input$grid_start_step2))
      
      n_durs2 <- length(seq(input$grid_dur_min2,
                            input$grid_dur_max2,
                            by = input$grid_dur_step2))
      
      n_total <- n_starts * n_durs + n_starts2 * n_durs2
      
    } else {
      
      n_total <- n_starts * n_durs
    }
    # n_total <- n_starts * n_durs
    
    HTML(sprintf(
      "<small><b>%d designs generated</b> (%d start offsets × %d durations) — reference: <b>%s</b></small>",
      n_total, n_starts, n_durs, input$grid_event
    ))
  })
  
  # ---- TAB 3:COMPUTATION----
  grid_results <- eventReactive(input$run_grid, {
    
    dt <- copy(dt_raw())
    req(dt)
    
    
    message(Sys.time(), "window tab - computation started")
    withProgress(message = "Building grid...", value = 0, {
      dt[, minute := hour(timestamp) * 60 + minute(timestamp)]
      dt[, date := as.Date(date)]
      
      comb <- unique(dt[, .(recorder, date)])
      
      start_offsets_min <- seq(input$grid_start_min * 60,
                               input$grid_start_max * 60,
                               by = input$grid_start_step * 60)
      
      durations_min <- seq(input$grid_dur_min * 60,
                           input$grid_dur_max * 60,
                           by = input$grid_dur_step * 60)
      
      solar_dt <- solar_table()
      setDT(solar_dt)
      message(Sys.time(), "window tab - starting calc richness")
      
      
      grid_all <- rbindlist(lapply(seq_len(nrow(comb)), function(i) {
        
        rec <- comb$recorder[i]
        d   <- comb$date[i]
        
        sun_row <- solar_dt[recorder == rec & date == d]
        if (nrow(sun_row) == 0) return(NULL)
        
        if (input$grid_event %in% c("dawn", "sunrise", "sunset", "dusk")) {
          ref <- sun_row[[input$grid_event]]
        }
        
        if (grepl("^\\d{2}:\\d{2}$", input$grid_event)) {
          ref <- as.POSIXct(paste(d, input$grid_event), tz = input$timezone)
        }
        
        ref_min <- hour(ref) * 60 + minute(ref)
        
        tmp <- as.data.table(expand.grid(
          recorder = rec,
          date = d,
          start_offset_min = start_offsets_min,
          duration = durations_min
        ))
        
        tmp[, start_abs_min := ref_min + start_offset_min]
        tmp
        
      }), use.names = TRUE)
      
      grid_all[, effort := duration]
      grid_all[, richness := NA_real_]
      setkey(dt, recorder, date)
      
      n <- nrow(grid_all)
      
      
      for (i in seq_len(n)) {
        
        day_data <- dt[.(grid_all$recorder[i], grid_all$date[i])]
        s0 <- grid_all$start_abs_min[i]
        s1 <- s0 + grid_all$duration[i]
        
        rec_i <- grid_all$recorder[i]
        d_i   <- grid_all$date[i]
        
        if (s0 >= 0 && s1 <= 1440) {
          
          # normal: within the day
          day_data <- dt[.(rec_i, d_i)]
          grid_all$richness[i] <- length(unique(day_data[minute >= s0 & minute <= s1, species]))
          
        } else if (s1 > 1440) {
          
          # crosses to next day: [s0..1440] on day d + [0..(s1-1440)] on day d+1
          day_data1 <- dt[.(rec_i, d_i)]
          day_data2 <- dt[.(rec_i, d_i + 1)]
          
          sp1 <- unique(day_data1[minute >= max(0, s0) & minute <= 1440, species])
          sp2 <- unique(day_data2[minute >= 0 & minute <= (s1 - 1440), species])
          
          grid_all$richness[i] <- length(unique(c(sp1, sp2)))
          
        } else if (s0 < 0) {
          
          # crosses to previous day: [ (1440+s0)..1440 ] on day d-1 + [0..s1] on day d
          day_data0 <- dt[.(rec_i, d_i - 1)]
          day_data1 <- dt[.(rec_i, d_i)]
          
          sp0 <- unique(day_data0[minute >= (1440 + s0) & minute <= 1440, species])
          sp1 <- unique(day_data1[minute >= 0 & minute <= min(1440, s1), species])
          
          grid_all$richness[i] <- length(unique(c(sp0, sp1)))
          
        } else {
          grid_all$richness[i] <- NA_real_
        }
        
        
      }
      message(Sys.time(), "window tab - big calc done ")
      # ---- build gris all----
      grid_all[, day_max := max(richness, na.rm = TRUE), by = date]
      grid_all[, richness_pct := 100 * richness / day_max]
      
      grid_all[, start_off_h := start_offset_min / 60]
      grid_all[, duration_h  := duration / 60]
      message(Sys.time(), "window tab - formatting result")
      grid_summary <- grid_all[, .(
        richness     = mean(richness, na.rm = TRUE),
        richness_pct = mean(richness_pct, na.rm = TRUE),
        effort       = mean(effort, na.rm = TRUE)
      ), by = .(start_off_h, duration_h)]
      message(Sys.time(), "window tab -2 ")
      # values$grid_results <- grid_summary
      grid_w1 <- copy(grid_summary)
      grid_w2 <- NULL
      grid_combined <- NULL
      message(Sys.time(), "window tab -3 ")
      if (isTRUE(input$use_window2_grid)) {
        message(Sys.time(), "window tab -4 ")
        start_offsets_min2 <- seq(input$grid_start_min2 * 60,
                                  input$grid_start_max2 * 60,
                                  by = input$grid_start_step2 * 60)
        message(Sys.time(), "window tab -5 ")
        durations_min2 <- seq(input$grid_dur_min2 * 60,
                              input$grid_dur_max2 * 60,
                              by = input$grid_dur_step2 * 60)
        
        # ---- build grid all - 2 window----
        grid_all2 <- rbindlist(lapply(seq_len(nrow(comb)), function(i) {
          
          rec <- comb$recorder[i]
          d   <- comb$date[i]
          
          sun_row <- solar_dt[recorder == rec & date == d]
          if (nrow(sun_row) == 0) return(NULL)
          
          if (input$grid_event2 %in% c("dawn", "sunrise", "sunset", "dusk")) {
            ref <- sun_row[[input$grid_event2]]
          }
          
          if (grepl("^\\d{2}:\\d{2}$", input$grid_event2)) {
            ref <- as.POSIXct(paste(d, input$grid_event2), tz = input$timezone)
          }
          
          ref_min <- hour(ref) * 60 + minute(ref)
          
          tmp <- as.data.table(expand.grid(
            recorder = rec,
            date = d,
            start_offset_min = start_offsets_min2,
            duration = durations_min2
          ))
          
          tmp[, start_abs_min := ref_min + start_offset_min]
          tmp
          
        }), use.names = TRUE)
        
        grid_all2[, effort := duration]
        grid_all2[, richness := NA_real_]
        
        # calc richness window 2 (même boucle que la tienne, mais sur grid_all2)
        n2 <- nrow(grid_all2)
        for (i in seq_len(n2)) {
          
          day_data <- dt[.(grid_all2$recorder[i], grid_all2$date[i])]
          
          s0 <- grid_all2$start_abs_min[i]
          s1 <- s0 + grid_all2$duration[i]
          
          grid_all2$richness[i] <- length(unique(
            day_data[minute >= s0 & minute <= s1, species]
          ))
        }
        
        grid_all2[, start_off_h := start_offset_min / 60]
        grid_all2[, duration_h  := duration / 60]
        
        grid_w2 <- grid_all2[, .(
          richness = mean(richness, na.rm = TRUE),
          effort   = mean(effort, na.rm = TRUE)
        ), by = .(start_off_h, duration_h)]
        # --- Combined: on pair les designs W1 et W2 par index (simple & stable)
        # base = grille W1, et on associe une ligne W2 via recyclage
        # (cohérent avec ton preview table cbind)
        design_w1 <- unique(grid_all[, .(start_abs_min, duration, start_off_h, duration_h)])
        design_w2 <- unique(grid_all2[, .(start_abs_min, duration, start_off_h, duration_h)])
        
        n1 <- nrow(design_w1)
        n2d <- nrow(design_w2)
        
        idx2 <- ((seq_len(n1) - 1) %% n2d) + 1L
        
        comb_design <- data.table(
          idx = seq_len(n1),
          start1 = design_w1$start_abs_min,
          dur1   = design_w1$duration,
          start_off_h = design_w1$start_off_h,
          duration_h  = design_w1$duration_h,
          start2 = design_w2$start_abs_min[idx2],
          dur2   = design_w2$duration[idx2]
        )
        message(Sys.time(), "window tab -5 ")
        # construire un grid_all_comb (recorder×date×design)
        grid_all_comb <- rbindlist(lapply(seq_len(nrow(comb)), function(i) {
          
          rec <- comb$recorder[i]
          d   <- comb$date[i]
          
          tmp <- copy(comb_design)
          tmp[, `:=`(recorder = rec, date = d)]
          tmp
        }))
        message(Sys.time(), "window tab -6 ")
        # grid_all_comb[, `:=`(effort = dur1 + dur2, richness = NA_real_)]
        grid_all_comb[, `:=`(
          effort = dur1 + dur2,
          richness = NA_real_,
          richness_w1 = NA_real_,
          richness_w2 = NA_real_
        )]
        message(Sys.time(), "window tab -7 ")
        dt_split <- split(dt, by = c("recorder", "date"), keep.by = FALSE)
        message(Sys.time(), "window tab -77 ")
        nC <- nrow(grid_all_comb)
        
        for (i in seq_len(nC)) {
          
          # day_data <- dt[.(grid_all_comb$recorder[i], grid_all_comb$date[i])]
          key <- paste(grid_all_comb$recorder[i], grid_all_comb$date[i], sep=".")
          day_data <- dt_split[[key]]
          if (is.null(day_data)) next
          s0a <- grid_all_comb$start1[i]
          s1a <- s0a + grid_all_comb$dur1[i]
          
          s0b <- grid_all_comb$start2[i]
          s1b <- s0b + grid_all_comb$dur2[i]
          sp1 <- unique(day_data[minute >= s0a & minute <= s1a, species])
          sp2 <- unique(day_data[minute >= s0b & minute <= s1b, species])
          
          grid_all_comb$richness_w1[i] <- length(sp1)
          grid_all_comb$richness_w2[i] <- length(sp2)
          grid_all_comb$richness[i]    <- length(unique(c(sp1, sp2)))
          
          if (i %% 100L == 0L || i == n) {
            incProgress(i / nC, detail = paste(i, "/", nC))
          }
          print(i)
        }
        message(Sys.time(), "window tab -777 ")
        grid_combined <- grid_all_comb[, .(
          richness = mean(richness, na.rm = TRUE),
          richness_w1 =  mean(richness_w1, na.rm = TRUE),
          richness_w2 =  mean(richness_w2, na.rm = TRUE),
          effort   = mean(effort, na.rm = TRUE)
        ), by = .(start_off_h, duration_h)]
      }
      
      # incProgress(0.5, detail = "half")
      message(Sys.time(), "window tab - ongoing")
      # --- Normalisation commune: 100% = max richesse sur la COMBINAISON
      if (isTRUE(input$use_window2_grid)) {
        message(Sys.time(), "window tab -8 ")
        max_global <- max(grid_combined$richness, na.rm = TRUE)
        
        grid_w1[, richness_pct := 100 * richness / max_global]
        grid_w2[, richness_pct := 100 * richness / max_global]
        grid_combined[, richness_pct := 100 * richness / max_global]
        grid_combined[, richness_w1_pct := 100 * richness_w1 / max_global]
        grid_combined[, richness_w2_pct := 100 * richness_w2 / max_global]
      } else {
        # comportement actuel (comme avant) — si tu veux garder ton "par date", laisse ton code.
        # mais pour être cohérent, on peut aussi faire:
        message(Sys.time(), "window tab -9")
        max_global <- max(grid_w1$richness, na.rm = TRUE)
        grid_w1[, richness_pct := 100 * richness / max_global]
      }
      message(Sys.time(), "window tab - results finalized")
      # incProgress(1, detail = "Done")
      values$grid_results <- list(w1 = grid_w1, w2 = grid_w2, comb = grid_combined)
      values$grid_results
      # grid_summary
    })
  })
  
  # ---- TAB 3: PLOT 1----
  
  output$grid_plot_w1 <- renderPlot({
    print(str(grid_results()))
    res <- grid_results()
    req(res$w1)
    
    df <- copy(res$w1)
    
    df[, richness_bin := floor(richness_pct / 5) * 5]
    frontier <- df[, .SD[which.min(duration_h)], by = richness_bin]
    frontier <- frontier[richness_bin > 0 & richness_bin < 100]
    message(Sys.time(), "window tab -9")
    dur_grid <- sort(unique(df$duration_h))
    
    p <- ggplot(df, aes(x = start_off_h, y = duration_h)) +
      geom_contour_filled(aes(z = richness_pct), color = "grey30", linewidth = 0.1, binwidth = 5) +
      scale_fill_scico_d(palette = "oslo", name = "Mean richness (%)") +
      geom_vline(
        xintercept = sort(unique(df$start_off_h)),
        color = "white", linewidth = 0.3, alpha = 0.5
      ) +
      geom_hline(
        yintercept = dur_grid,
        color = "white", linewidth = 0.3, alpha = 0.5
      ) +
      labs(
        x = paste0("Start offset relative to ", input$grid_event, " (hours)"),
        y = "Recording duration (hours)", title = "Window 1 "
      ) +
      theme_minimal() +
      theme(panel.grid = element_blank())
    
    if (isTRUE(input$show_frontier)) {
      p <- p + geom_path(
        data = frontier[order(richness_bin)],
        aes(x = start_off_h, y = duration_h),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.9
      )
    }
    
    p
  })
  
  # ---- TAB 3: PLOT 2----
  output$grid_plot_w2 <- renderPlot({
    req(input$use_window2_grid)
    res <- grid_results()
    req(res$w2)
    message(Sys.time(), "window tab -10 ")
    df <- copy(res$w2)
    df[, richness_bin := floor(richness_pct / 5) * 5]
    frontier <- df[, .SD[which.min(duration_h)], by = richness_bin]
    frontier <- frontier[richness_bin > 0 & richness_bin < 100]
    
    dur_grid <- sort(unique(df$duration_h))
    
    p <- ggplot(df, aes(x = start_off_h, y = duration_h)) +
      geom_contour_filled(aes(z = richness_pct), color = "grey30", linewidth = 0.1, binwidth = 5) +
      scale_fill_scico_d(palette = "bilbao", name = "Mean richness (%)") +  # <- couleur différente
      geom_vline(xintercept = sort(unique(df$start_off_h)), color = "white", linewidth = 0.3, alpha = 0.5) +
      geom_hline(yintercept = dur_grid, color = "white", linewidth = 0.3, alpha = 0.5) +
      labs(
        x = paste0("Start offset relative to ", input$grid_event2, " (hours)"),
        y = "Recording duration (hours)", title = "Window 2 "
      ) +
      theme_minimal() +
      theme(panel.grid = element_blank())
    
    if (isTRUE(input$show_frontier)) {
      p <- p + geom_path(
        data = frontier[order(richness_bin)],
        aes(x = start_off_h, y = duration_h),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.9
      )
    }
    
    p
  })
  
  # ---- TAB 3: PLOT COMB----
  output$grid_plot_combined <- renderPlot({
    
    res <- grid_results()
    req(res$w1, res$w2, res$comb)
    
    
    df1 <- copy(res$w1)
    df1[, type := "Window 1"]
    df1[, effort_total := duration_h]
    
    df2 <- copy(res$w2)
    df2[, type := "Window 2"]
    df2[, effort_total := duration_h]
    
    df3 <- copy(res$comb)
    df3[, type := "Combined"]
    df3[, effort_total := effort / 60]   # effort était en minutes
    
    df <- rbindlist(list(df1, df2, df3), use.names = TRUE, fill = TRUE)
    # --- segments reliant W1 et W2 au Combined ---
    seg_w1 <- df3[, .(
      x  = duration_h,
      y  = richness_w1_pct,
      xend = effort_total,
      yend = richness_pct
    )]
    
    seg_w2 <- df3[, .(
      x  = duration_h,
      y  = richness_w2_pct,
      xend = effort_total,
      yend = richness_pct
    )]
    
    ggplot() +
      
      # segments W1 → Combined
      geom_segment(
        data = seg_w1,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = "grey60", alpha = 0.8, linewidth = 0.3
      ) +
      geom_segment(
        data = seg_w2,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = "grey60", alpha = 0.8, linewidth = 0.3
      ) +
      
      # points Window 1
      geom_point(
        data = df1,
        aes(x = effort_total, y = richness_pct),
        color = "steelblue",
        size = 3,
        alpha = 0.9
      ) +
      
      # points Window 2
      geom_point(
        data = df2,
        aes(x = effort_total, y = richness_pct),
        color = "orange",
        size = 3,
        alpha = 0.9
      ) +
      
      # points Combined (plus gros)
      geom_point(
        data = df3,
        aes(x = effort_total, y = richness_pct),
        color = "darkgreen",
        size = 4.5
      ) +
      
      labs(
        x = "Total recording effort (hours)",
        y = "Mean relative richness (%)"
      ) +
      theme_minimal()
  })
  # ---- plot organisation----
  output$grid_plots_ui <- renderUI({
    if (isTRUE(input$use_window2_grid)) {
      tagList(
        plotOutput("grid_plot_w1", height = "450px"),
        plotOutput("grid_plot_w2", height = "450px"),
        plotOutput("grid_plot_combined", height = "450px")
      )
    } else {
      plotOutput("grid_plot_w1", height = "650px")
    }
  })
  
  # ---- TAB 4: MULTI OPTI----
  multi_combo <- eventReactive(input$run_multi, {
    
    req(values$summary_richness)
    req(values$grid_results)
    req(values$schedules)
    message(Sys.time(), "multi optimum- starting")
    showNotification("starting analysis", type = "message", duration = 2)
    
    summary_richness <- values$summary_richness
    # grid_results     <- values$grid_results
    schedules        <- values$schedules
    grid_results_raw <- values$grid_results
    print(names(grid_results_raw))
    print(str(grid_results_raw))
    req(!is.null(grid_results_raw))
    grid_results_raw <- values$grid_results
    req(!is.null(grid_results_raw))
    
    if (is.list(grid_results_raw)) {
      
      two_windows <- !is.null(grid_results_raw$w2)
      
    } else {
      
      two_windows <- FALSE
    }
    
    
    if (two_windows) {
      
      window_sets <- list(
        W1       = grid_results_raw$w1,
        W2       = grid_results_raw$w2,
        Combined = grid_results_raw$comb
      )
      
    } else {
      
      window_sets <- list(
        W1 = if (is.list(grid_results_raw)) grid_results_raw$w1 else grid_results_raw
      )
    }
    
    message(Sys.time(), "multi optimum- 3")
    dt <- copy(dt_raw())
    req(dt)
    dt[, date := as.Date(date)]
    solar_dt <- solar_table()
    setDT(solar_dt)
    message(Sys.time(), "multi optimum- starting loop")
    # print(class(grid_results)); print(nrow(grid_results))
    withProgress(message = "Running multi-optimum analysis...", value = 0, {
      
      combo_results <- list()
      k <- 1
      message(Sys.time(), "multi optimum- 4")
      # ---- TOTAL NUMBER OF ITERATIONS ----
      total_iter <- 0L
      
      for (sn in names(window_sets)) {
        gr <- window_sets[[sn]]
        if (!is.null(gr)) {
          total_iter <- total_iter + nrow(gr)
        }
      }
      
      done_iter <- 0L
      for (set_name in names(window_sets)) {
        
        grid_results <- window_sets[[set_name]]
        message(Sys.time(), "multi optimum- 44")
        
        if (is.null(grid_results)) next
        if (!is.data.table(grid_results)) setDT(grid_results)
        nW <- nrow(grid_results)
        if (nW == 0) next
        
        for (w in seq_len(nW)) {
          
          message(Sys.time(), "multi optimum- 5")
          done_iter <- done_iter + 1L
          
          if (done_iter %% 10L == 0L || done_iter == total_iter) {
            incProgress(done_iter / total_iter,
                        detail = paste(done_iter, "/", total_iter))
          }
          
          
          if (set_name == "Combined") {
            win_row_w1 <- window_sets$W1[w]
            win_row_w2 <- window_sets$W2[w]
            
            w1_start_off <- win_row_w1$start_off_h
            w1_dur_h     <- win_row_w1$duration_h
            
            w2_start_off <- win_row_w2$start_off_h
            w2_dur_h     <- win_row_w2$duration_h
            
          } else if (set_name == "W1") {
            win_row <- grid_results[w]
            w1_start_off <- win_row$start_off_h
            w1_dur_h     <- win_row$duration_h
            w2_start_off <- NA_real_
            w2_dur_h     <- NA_real_
            
          } else if (set_name == "W2") {
            win_row <- grid_results[w]
            w1_start_off <- NA_real_
            w1_dur_h     <- NA_real_
            w2_start_off <- win_row$start_off_h
            w2_dur_h     <- win_row$duration_h
          }
          
          message(Sys.time(), "multi optimum- 6")
          
          for (p in seq_along(schedules)) {
            
            message(Sys.time(), "multi optimum- 7")
            duty_sched <- schedules[[p]]
            sched_combined_list <- list()
            
            comb_rd <- unique(duty_sched[, .(recorder, date)])
            
            for (j in seq_len(nrow(comb_rd))) {
              
              message(Sys.time(), "multi optimum- 8")
              rec <- comb_rd$recorder[j]
              d   <- comb_rd$date[j]
              
              sun_row <- solar_dt[recorder == rec & date == d]
              if (nrow(sun_row) == 0) next
              message(Sys.time(), "multi optimum- 88")
              
              # ---- ref_min window  1 ----
              if (input$grid_event %in% c("dawn","sunrise","sunset","dusk")) {
                ref1 <- sun_row[[input$grid_event]]
              } else if (grepl("^\\d{2}:\\d{2}$", input$grid_event)) {
                ref1 <- as.POSIXct(paste(d, input$grid_event), tz = input$timezone)
              }
              ref1_min <- hour(ref1) * 60 + minute(ref1)
              
              # ---- ref_min window 2 ----
              # 
              if (!is.null(input$grid_event2)) {
                if (input$grid_event2 %in% c("dawn","sunrise","sunset","dusk")) {
                  ref2 <- sun_row[[input$grid_event2]]
                } else if (grepl("^\\d{2}:\\d{2}$", input$grid_event2)) {
                  ref2 <- as.POSIXct(paste(d, input$grid_event2), tz = input$timezone)
                }
                ref2_min <- hour(ref2) * 60 + minute(ref2)
              } else {
                # fallback si pas de window2 event défini (évite que tout devienne identique)
                ref2_min <- ref1_min
              }
              
              # ---- construire fenêtres ----
              if (!is.na(w1_start_off)) {
                w1_start <- ref1_min + w1_start_off * 60
                w1_end   <- w1_start + w1_dur_h * 60
              }
              if (!is.na(w2_start_off)) {
                w2_start <- ref2_min + w2_start_off * 60
                w2_end   <- w2_start + w2_dur_h * 60
              }
              
              sched_day <- duty_sched[recorder == rec & date == d]
              message(Sys.time(), "multi optimum- 888")
              
              if (set_name == "W1") {
                message(Sys.time(), "multi optimum- 9")
                sched_day <- sched_day[start_min >= w1_start & end_min <= w1_end]
                message(Sys.time(), "multi optimum- 99")
                
              } else if (set_name == "W2") {
                message(Sys.time(), "multi optimum- 999")
                sched_day <- sched_day[start_min >= w2_start & end_min <= w2_end]
                message(Sys.time(), "multi optimum- 9999")
                
              } else if (set_name == "Combined") {
                message(Sys.time(), "multi optimum- 99a")
                sched_day <- sched_day[
                  (start_min >= w1_start & end_min <= w1_end) |
                    (start_min >= w2_start & end_min <= w2_end)
                ]
                message(Sys.time(), "multi optimum- 99b")
              }
              
              sched_combined_list[[length(sched_combined_list) + 1]] <- sched_day
            }
            
            sched_combined <- rbindlist(sched_combined_list)
            message(Sys.time(), "multi optimum- 85")
            if (nrow(sched_combined) == 0) next
            
            message(Sys.time(), "multi optimum- 9")
            res <- compute_richness(
              dt = dt,
              schedule = sched_combined,
              per_spatial = "all",
              per_temporal = "day"
            )
            message(Sys.time(), "multi optimum- 10")
            
            combo_results[[k]] <- data.table(
              window_type = set_name,
              start1_hour = if (set_name %in% c("W1","Combined")) w1_start_off else NA_real_,
              dur1_h      = if (set_name %in% c("W1","Combined")) w1_dur_h else NA_real_,
              start2_hour = if (set_name %in% c("W2","Combined")) w2_start_off else NA_real_,
              dur2_h      = if (set_name %in% c("W2","Combined")) w2_dur_h else NA_real_,
              duty_period = attr(duty_sched, "period"),
              effort      = sum(sched_combined$end_min - sched_combined$start_min + 1),
              richness    = ifelse(length(res$richness) == 0, 0, mean(res$richness, na.rm = TRUE))
            )
            
            k <- k + 1
          }
        }
      }
      message(Sys.time(), "multi optimum - loop done")
      print(unique(combo_results[[1]]$window_type))
      combo <- rbindlist(combo_results)
      combo[, richness_pct := 100 * richness / max(richness)]
      print(table(combo$window_type))
      print(combo[, .(mean_eff = mean(effort),
                      mean_rich = mean(richness)), 
                  by = window_type])
      if (nrow(combo) == 0 || all(is.na(combo$richness))) {
        warning("Combo empty — no valid schedules")
        return(NULL)
      }
      
      combo[, pareto := FALSE]
      best <- -Inf
      for (i in seq_len(nrow(combo))) {
        if (!is.na(combo$richness_pct[i]) && combo$richness_pct[i] > best) {
          combo$pareto[i] <- TRUE
          best <- combo$richness_pct[i]
        }
      }
      message(Sys.time(), "multi optimum- caculation done result")
      combo
    })
  })
  
  #---- TAB 4 PLOT ----
  output$plot_multi_optimum <- renderPlotly({
    
    dt <- multi_combo()
    req(dt)
    message(Sys.time(), "multi optimum- starting plot")
    dt[, dominated := FALSE]
    
    for (i in seq_len(nrow(dt))) {
      dominated_flag <- any(
        dt$effort <= dt$effort[i] &
          dt$richness_pct >= dt$richness_pct[i] &
          (dt$effort < dt$effort[i] |
             dt$richness_pct > dt$richness_pct[i])
      )
      dt$dominated[i] <- dominated_flag
    }
    message(Sys.time(), "multi optimum- pareto calc done")
    dt[, pareto := !dominated]
    pareto_front <- dt[pareto == TRUE][order(effort)]
    
    durations <- battery_specs$duration_h
    codes     <- battery_specs$code
    
    message(Sys.time(), "multi optimum- format done, start plot")
    
    p <- ggplot(dt, aes(x = effort, y = richness_pct)) +
      geom_point(
        aes(
          color = window_type,
          shape = window_type,
          
          text = paste0(
            "Richness: ", round(richness_pct, 1), "%<br>",
            "<b>Type:</b> ", window_type, "<br>",
            "Effort: ", effort, " min - ",  round(.data$effort / 60, 2), " hours ", "<br>",
            ifelse(!is.na(start2_hour),
                   paste0(
                     "<b>Window 1</b><br>",
                     "Start offset: ", round(start1_hour, 2), " h<br>",
                     "Duration: ", round(dur1_h, 2), " h<br>",
                     "<b>Window 2</b><br>",
                     "Start offset: ", round(start2_hour, 2), " h<br>",
                     "Duration: ", round(dur2_h, 2), " h<br>"
                   ),
                   paste0(
                     "<b>Window</b><br>",
                     "Start offset: ", round(start1_hour, 2), " h<br>",
                     "Duration: ", round(dur1_h, 2), " h<br>"
                   )
            ),
            # "Start: ", start_hour_1, "<br>",
            # "Duration: ", dur1_h, "<br>",
            "Period: ", duty_period
          )
        ),
        alpha = 0.7,
        size = 2
      ) +
      geom_line(
        data = pareto_front,
        aes(group = 1),
        color = "red",
        linewidth = 1
      ) +
      geom_point(
        data = pareto_front,
        aes(text = paste0(
          "Local optimum<br>",
          "Richness: ", round(richness_pct, 1), "%<br>",
          "<b>Type:</b> ", window_type, "<br>",
          "Effort: ", effort, " min - ",  round(.data$effort / 60, 2), " hours ", "<br>",
          ifelse(!is.na(start2_hour),
                 paste0(
                   "<b>Window 1</b><br>",
                   "Start offset: ", round(start1_hour, 2), " h<br>",
                   "Duration: ", round(dur1_h, 2), " h<br>",
                   "<b>Window 2</b><br>",
                   "Start offset: ", round(start2_hour, 2), " h<br>",
                   "Duration: ", round(dur2_h, 2), " h<br>"
                 ),
                 paste0(
                   "<b>Window</b><br>",
                   "Start offset: ", round(start1_hour, 2), " h<br>",
                   "Duration: ", round(dur1_h, 2), " h<br>"
                 )
          ),
          # "Start: ", start_hour, "<br>",
          # "Duration: ", duration_h, "<br>",
          "Period: ", duty_period,
          "<br><br><b>Full battery sets needed:</b><br>",
          paste0(
            codes[1], ": ", round((.data$effort / 60) / durations[1], 2), "<br>",
            codes[2], ": ", round((.data$effort / 60) / durations[2], 2), "<br>",
            codes[3], ": ", round((.data$effort / 60) / durations[3], 2), "<br>",
            codes[4], ": ", round((.data$effort / 60) / durations[4], 2), "<br>",
            codes[5], ": ", round((.data$effort / 60) / durations[5], 2), "<br>",
            codes[6], ": ", round((.data$effort / 60) / durations[6], 2)
          )
        )),
        color = "red",
        size = 2
      ) +
      scale_color_manual(
        values = c(
          W1 = "#1f78b4",
          W2 = "#6a3d9a",
          Combined = "#33a02c"
        ),
        name = "Window type"
      ) +
      scale_shape_manual(
        values = c(
          W1 = 16,
          W2 = 17,
          Combined = 18
        ),
        name = "Window type"
      )
    # scale_color_viridis_c() +
    labs(
      x = "Recording effort (min)",
      y = "Relative richness (%)",
      title = "Pareto frontier: optimal trade-off richness vs effort"
    ) +
      theme_minimal()
    
    ggplotly(p, tooltip = "text")
  })
  
  # ---- TAB 5: TARGET SPECIES----
  
  observe({
    req(dt_raw())
    message(Sys.time(), "target species - choices")
    updateSelectizeInput(
      session,
      "target_species",
      choices = sort(unique(dt_raw()$species)),
      server = TRUE
    )
  })
  
  target_result <- eventReactive(input$run_target, {
    
    req(dt_raw())
    req(input$target_species)
    message(Sys.time(), "target species - starting")
    dt <- copy(dt_raw())
    recs <- unique(dt$recorder)
    
    solar_dt <- solar_table()
    solar_dt <- setDT(solar_dt)
    
    sp <- dt[species %in% input$target_species]
    sp[, date := as.Date(date)]
    message(Sys.time(), "target species - parameter")
    time_bin_min <- as.numeric(input$time_resolution)
    sp[, time_bin := (minute_of_day %/% time_bin_min) * time_bin_min]
    
    time_levels <- seq(0, 24 * 60 - time_bin_min, by = time_bin_min)
    message(Sys.time(), "target species - building grid")
    grid <- CJ(
      date = seq(min(sp$date), max(sp$date), by = "day"),
      time_bin = time_levels
    )
    
    periods <- seq(
      input$target_period_min,
      input$target_period_max,
      by = input$target_period_step
    )
    message(Sys.time(), "target species - precomputation done, start loop")
    all_heatmaps <- list()
    
    withProgress(message = "Computing activity...", value = 0, {
      
      for (i in seq_along(periods)) {
        
        p <- periods[i]
        incProgress(1 / length(periods), detail = paste(i, "/", length(periods)))
        
        sched_rec_list <- lapply(recs, function(r) {
          
          sun_r <- solar_dt[solar_dt$recorder == r, ]
          
          sched_r <- create_sun_duty_schedule(
            sun_table = sun_r,
            period = p,
            duty_duration = input$target_duty_duration,
            nb_duty = input$target_nb_duty,
            start_event = "dawn",
            end_event   = "dusk",
            extend_before_hours = 1,
            extend_after_hours  = 1
          )
          sched_r[, recorder := r]
          sched_r
        })
        
        sched <- rbindlist(sched_rec_list)
        
        sp_sched <- sched[
          sp,
          on = .(date,
                 start_min <= minute_of_day,
                 end_min >= minute_of_day),
          nomatch = 0
        ]
        
        counts <- sp_sched[, .(n = .N), by = .(date, time_bin)]
        
        counts <- merge(grid, counts, by = c("date", "time_bin"), all.x = TRUE)
        counts[is.na(n), n := 0]
        
        counts[, period_label := paste0("1/", p)]
        counts[, period_num := p]
        
        all_heatmaps[[as.character(p)]] <- counts
      }
    })
    message(Sys.time(), "target species - loop done ")
    heatmap_dt <- rbindlist(all_heatmaps)
    
    heatmap_dt[, period_label :=
                 factor(period_label,
                        levels = paste0("1/", sort(unique(period_num))))]
    
    heatmap_dt[, date_chr := as.character(date)]
    
    heatmap_dt[, hover_text := paste0(
      "Time: ",
      sprintf("%02d:%02d",
              time_bin %/% 60,
              time_bin %% 60),
      "<br>Date: ", date_chr,
      "<br>Activity: ", n, " detections"
    )]
    message(Sys.time(), "target species - resukt finished ")
    heatmap_dt
  })
  
  # ---- TAB 5: PLOT ----
  output$target_heatmap <- renderPlotly({
    
    heatmap_dt <- target_result()
    req(heatmap_dt)
    message(Sys.time(), "starting plot - ")
    global_max <- max(heatmap_dt$n, na.rm = TRUE)
    
    pal_main <- hcl.colors(1000, "Berlin")
    zero_col <- "grey89"
    pal <- c(zero_col, pal_main)
    
    fig <- plot_ly(
      data = heatmap_dt,
      x = ~date_chr,
      y = ~time_bin / 60,
      z = ~n,
      type = "heatmap",
      frame = ~period_label,
      colors = pal,
      zmin = 0,
      zmax = global_max,
      text = ~hover_text,
      hovertemplate = "%{text}<extra></extra>",
      xgap = 0.3,
      ygap = 0.3,
      showscale = TRUE,
      opacity = 0.9
    )
    
    fig %>%
      layout(
        title = paste("Species detection for ", input$target_species),
        yaxis = list(
          title = "Hour of day",
          tickmode = "linear",
          tick0 = 0,
          dtick = 1,
          range = c(0, 24),
          showgrid = FALSE
        ),
        xaxis = list(
          title = "Date",
          showgrid = FALSE,
          tickfont = list(size = 10),
          automargin = TRUE
        )
      ) %>%
      config(
        displaylogo = FALSE,
        modeBarButtonsToRemove = c("select2d", "lasso2d")
      )
  })
}


# ---- Launch ----
shinyApp(ui, server)
