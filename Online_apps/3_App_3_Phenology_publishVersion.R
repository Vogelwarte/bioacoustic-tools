# ==============================================================================
# BirdNET-ResChecker - PHENOLOGY (online version)
# ------------------------------------------------------------------------------
# Loads BirdNET results (selection tables or a compiled file), converts the
# recorder timestamps to the correct time zone and draws a phenology heatmap
# (detections per time slot and per day, with sunrise / sunset / dawn / dusk).
# ==============================================================================


# --- 1. PACKAGES --------------------------------------------------------------
library(shiny)
library(bslib)
library(stringr)
library(stringi)
library(fs)
library(data.table)
library(ggplot2)
library(plotly)
library(hms)
library(suncalc)
library(scales)
library(lubridate)
# dipsaus is no longer needed: the folder upload uses a standard fileInput

# Allow large uploads (default Shiny limit is 5 MB). Adjust to the server RAM.
options(shiny.maxRequestSize = 500 * 1024^2)


# --- 2. CUSTOM FUNCTIONS ------------------------------------------------------
source("Common_functions/pheno_matrix_plotly.R")
source("Common_functions/clean_text.R")
source("Common_functions/load_selection_tables_ONLINE.R")


# --- 3. HELPER: TIME ZONES, DATE STRINGS AND RECORDER NAMES -------------------
# Takes the table returned by the loader and:
#   1. creates 'date_strings' (YYYYMMDD_HHMMSS) from the audio file name if missing
#   2. creates 'recorder' from the file name if missing
#   3. reads the recorder clock in 'device_tz' (gives the correct instant)
#   4. expresses that same instant in 'deployment_tz' (local time)
#   5. adds the offset of each detection within the file (Start/Stop_segment)
#   6. derives Date, Hour, Min and Time in local time
apply_timezones <- function(DT, device_tz, deployment_tz, compiled) {
  DT <- data.table::as.data.table(data.table::copy(DT))
  
  # Columns that may contain the audio file name, depending on the format
  file_cols <- intersect(
    c("Begin_Path", "Begin.Path", "Begin.File", "File", "file", "Source_file"),
    names(DT)
  )
  
  # 1. date_strings: ALWAYS read again from the audio file name
  #    (the loader may have built them in the server's time zone)
  if (length(file_cols) > 0) {
    fname <- basename(as.character(DT[[file_cols[1]]]))
    ds <- stringr::str_extract(fname, "\\d{8}_\\d{6}")
    if (!all(is.na(ds))) DT[, date_strings := ds]
  }
  if (!"date_strings" %in% names(DT) || all(is.na(DT$date_strings))) {
    stop("Could not read a date (YYYYMMDD_HHMMSS) from the file names. Columns are: ",
         paste(names(DT), collapse = ", "))
  }
  
  # 2. recorder = part of the file name before the date (e.g. "REC01")
  if (!"recorder" %in% names(DT) && length(file_cols) > 0) {
    fname <- basename(as.character(DT[[file_cols[1]]]))
    rec <- sub("_?\\d{8}_\\d{6}.*$", "", fname)
    rec[rec == "" | rec == fname] <- NA
    DT[, recorder := rec]
  }
  
  # 3. Recorder clock -> correct instant
  DT[, DateTime_Display := as.POSIXct(date_strings, format = "%Y%m%d_%H%M%S", tz = device_tz)]
  
  # 4. Same instant, shown in local (deployment) time
  DT[, DateTime_Real := lubridate::with_tz(DateTime_Display, tzone = deployment_tz)]
  
  # 5. Offset of the detection within the audio file (seconds)
  if (compiled && "File.Offset..s." %in% names(DT)) {
    DT[, time_offset := File.Offset..s.]
  } else if ("Begin.Time..s." %in% names(DT)) {
    DT[, time_offset := Begin.Time..s.]
  } else {
    DT[, time_offset := 0]
  }
  
  DT[, Start_segment := DateTime_Real + time_offset]
  
  if ("End.Time..s." %in% names(DT)) {
    if (compiled && "File.Offset..s." %in% names(DT)) {
      DT[, Stop_segment := Start_segment + (End.Time..s. - File.Offset..s.)]
    } else {
      DT[, Stop_segment := DateTime_Real + End.Time..s.]
    }
  } else {
    DT[, Stop_segment := Start_segment]
  }
  
  # 6. Derived columns in local time (tz is required, otherwise as.Date uses UTC)
  DT[, `:=`(
    Date = as.Date(Start_segment, tz = deployment_tz),
    Hour = lubridate::hour(Start_segment),
    Min  = lubridate::minute(Start_segment),
    Time = sprintf("%02d:%02d", lubridate::hour(Start_segment), lubridate::minute(Start_segment))
  )]
  DT[]
}


# --- 4. GLOBAL SETTINGS -------------------------------------------------------
# Time zone choices: fixed offsets first (for recorders that ignore DST),
# then all named time zones.
fixed_tz <- c(
  "UTC"           = "UTC",
  "UTC+1 (fixed)" = "Etc/GMT-1",
  "UTC+2 (fixed)" = "Etc/GMT-2",
  "UTC+3 (fixed)" = "Etc/GMT-3",
  "UTC-1 (fixed)" = "Etc/GMT+1",
  "UTC-2 (fixed)" = "Etc/GMT+2"
)
all_tz     <- OlsonNames()
tz_choices <- c(fixed_tz, setNames(all_tz, all_tz))

# Dark theme
dark_theme <- bs_theme(
  version   = 5,
  bootswatch = "darkly",
  bg        = "#1a1a1a",
  fg        = "#e0e0e0",
  primary   = "#375a7f",
  secondary = "#444444"
)


# --- 5. USER INTERFACE --------------------------------------------------------
ui <- page_sidebar(
  theme = dark_theme,
  
  # ---- Sidebar: data source, time zones, recorders ----
  sidebar = sidebar(
    title = "Data Source",
    
    checkboxInput("Compiled_F", "Results in compiled format?", value = FALSE),
    
    # Folder upload (selection tables)
    # Standard Shiny fileInput with the 'webkitdirectory' attribute: the browser
    # opens a FOLDER picker and Shiny uploads every file of that folder, with
    # its own progress bar. Works in Chrome, Edge, Firefox and Safari.
    conditionalPanel(
      condition = "input.Compiled_F == false",
      tagAppendAttributes(
        fileInput("dir1", "Choose results folder",
                  multiple = TRUE, buttonLabel = "Browse folder..."),
        webkitdirectory = NA, directory = NA,
        .cssSelector = "#dir1"   # the <input type="file"> carries the input id
      ),
      helpText("Pick the folder with the BirdNET result tables. Only .txt / .csv ",
               "files are uploaded (audio is skipped). Confirm the browser prompt ",
               "and wait for 'Upload complete'.")
    ),
    
    # Single-file upload (compiled results)
    conditionalPanel(
      condition = "input.Compiled_F == true",
      fileInput("compiled_file", "Choose compiled results file",
                multiple = FALSE, accept = c(".csv", ".txt"))
    ),
    
    actionButton("start", "Load data", class = "btn-primary", width = "100%"),
    hr(),
    
    selectInput("device_tz", "Time zone of the recorder clock",
                choices = tz_choices, selected = "Etc/GMT-1"),
    selectInput("deployment_tz", "Time zone of the deployment site (local time)",
                choices = tz_choices, selected = "Europe/Zurich"),
    helpText("Set the time zones BEFORE loading the data."),
    checkboxInput("fixed_recorder_clock",
                  "Plot on recorder clock (no daylight saving time)", value = TRUE),
    actionButton("restart_timezone", "Update time zones", class = "btn-primary", width = "100%"),
    hr(),
    
    uiOutput("recorder_ui"),
    hr(),
    
    textOutput("status")
  ),
  
  # ---- Header with logo ----
  div(
    style = "display: flex; justify-content: space-between; align-items: center; padding: 10px 20px; background-color: #0f0f0f; border-bottom: 1px solid #333; margin-bottom: 20px;",
    div(style = "font-size: 24px; font-weight: bold; color: #fff;",
        "BirdNET-ResChecker - PHENOLOGY"),
    div(tags$img(src = "logo.png", height = "70px",
                 style = "margin-right: 10px; border-radius: 12px; padding: 6px 10px; background-color: white;"))
  ),
  
  # ---- Main panels ----
  navset_card_underline(
    nav_panel(
      "Phenology",
      layout_columns(
        col_widths = c(4, 8),
        
        # Plot parameters
        card(
          full_screen = TRUE,
          card_header("Phenology parameters"),
          selectizeInput("species_pheno", "Species",
                         choices = NULL, selected = "All species", multiple = TRUE,
                         options = list(placeholder = 'Select species or "All species"', create = FALSE)),
          numericInput("Confid_Pheno", "Confidence threshold", min = 0.01, max = 1, value = 0.01),
          layout_columns(
            col_widths = c(6, 6),
            numericInput("Lat", "Latitude",  min = -90,  max = 90,  value = 46.97),
            numericInput("Lon", "Longitude", min = -180, max = 180, value = 6.97)
          ),
          sliderInput("Unit", "Aggregation interval (min)", min = 1, max = 60, value = 15),
          checkboxInput("Noctu_plot", "Nocturnal plot (noon to noon)", value = FALSE),
          checkboxInput("log_colour", "Log colour scale (shows rare and frequent slots)", value = TRUE),
          sliderInput("tile_alpha", "Heatmap opacity",
                      min = 0.1, max = 1, value = 0.9, step = 0.05, width = "100%"),
          actionButton("start_pheno", "Generate plot", class = "btn-primary", width = "100%")
        ),
        
        # Plot
        card(
          full_screen = TRUE,
          card_header("Phenology plot"),
          plotlyOutput("pheno_plot", height = "700px")
        )
      )
    ),
    
    nav_panel(
      "Reference",
      layout_columns(
        col_widths = c(2, 8, 2),
        textOutput("description"),
        textOutput("contributions"),
        textOutput("license")
      )
    )
  ),
  
  # ---- Custom CSS for the dark theme ----
  tags$head(
    # Folder upload in ONE request.
    # Shiny's fileInput uploads files one by one (one round trip per file),
    # which is very slow for folders with many small tables. This script
    # catches the folder selection BEFORE Shiny sees it (capture phase), packs
    # the .txt / .csv tables into a single zip in the browser (JSZip) and
    # hands only that zip to Shiny. The server unzips it (section 6.3).
    # If JSZip cannot be loaded, the normal file-by-file upload is used.
    # Offline server? Save jszip.min.js in www/ and use src = "jszip.min.js".
    tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/jszip/3.10.1/jszip.min.js"),
    tags$script(HTML("
      document.addEventListener('change', function (e) {
        var el = e.target;
        if (!el || el.id !== 'dir1' || !el.files || el.files.length === 0) return;
        if (el._bundled) { el._bundled = false; return; }   // our zip: let Shiny upload it

        // keep only the result tables
        var files = Array.prototype.filter.call(el.files, function (f) {
          return /\\.(txt|csv)$/i.test(f.name);
        });
        if (typeof JSZip === 'undefined') {               // fallback: filter only
          var dt0 = new DataTransfer();
          files.forEach(function (f) { dt0.items.add(f); });
          el.files = dt0.files;
          return;
        }
        e.stopImmediatePropagation();                     // hold back Shiny's handler
        if (files.length === 0) { console.log('No .txt / .csv files in this folder'); return; }

        var box = $(el).closest('.input-group').find('input[type=text]');
        box.val('Packing ' + files.length + ' files...');

        var zip = new JSZip(), used = {};
        files.forEach(function (f) {
          var nm = f.name, k = 1;
          while (used[nm]) nm = (k++) + '_' + f.name;     // same name in two sub-folders
          used[nm] = true;
          zip.file(nm, f);
        });
        zip.generateAsync({ type: 'blob', compression: 'DEFLATE',
                            compressionOptions: { level: 3 } })
          .then(function (blob) {
            var dt = new DataTransfer();
            dt.items.add(new File([blob], 'results_folder.zip', { type: 'application/zip' }));
            el.files = dt.files;
            el._bundled = true;
            el.dispatchEvent(new Event('change', { bubbles: true }));  // now Shiny uploads
          })
          .catch(function (err) { box.val('Packing failed: ' + err); });
      }, true);
    ")),
    tags$style(HTML("
      body { background-color: #2b2b2b; color: #e0e0e0; }
      .form-control, .selectize-input, .selectize-control.multi .selectize-input > div {
        background-color: #444444 !important; color: #ffffff !important; border: 1px solid #555555 !important;
      }
      .selectize-dropdown, .selectize-dropdown-content { background-color: #444444 !important; color: #ffffff !important; }
      .selectize-dropdown .option { color: #ffffff !important; }
      .selectize-dropdown .active { background-color: #555555 !important; color: #ffffff !important; }
      .control-label, .bslib-sidebar-input label { color: #ffffff !important; font-weight: 600; }
      .card { background-color: #333333; border: 1px solid #444444; color: #ffffff; }
      .card-header { background-color: #3a3a3a; border-bottom: 1px solid #444444; color: #ffffff; }
      ::-webkit-scrollbar { width: 10px; }
      ::-webkit-scrollbar-track { background: #2b2b2b; }
      ::-webkit-scrollbar-thumb { background: #555555; border-radius: 5px; }
      ::-webkit-scrollbar-thumb:hover { background: #777777; }
    "))
  )
)


# --- 6. SERVER ----------------------------------------------------------------
server <- function(input, output, session) {
  
  # ---- 6.1 Static texts (Reference tab) ----
  output$description   <- renderText("This application analyses BirdNET detection results and produces phenology plots. Upload your results on the left, load the data, then generate the plot.")
  output$contributions <- renderText("Contributions: Amandine Serrurier, Jean-Nicolas Pradervand & Christophe Sahli\nSwiss Ornithological Institute")
  output$license       <- renderText("MIT License \u00a9 2025 Christophe Sahli")  
  # ---- 6.2 Reactive storage ----
  DT_reac <- reactiveVal(NULL)   # full, time-corrected detection table
  dir1    <- reactiveVal(NULL)   # server folder holding the uploaded files
  
  output$status <- renderText("No data loaded.")
  
  # ---- 6.3 Folder upload: rebuild the folder on the server ----
  # Shiny stores uploads as 0.txt, 1.txt, ... in separate temp folders.
  # The loader expects one folder with the original file names, so the files
  # are copied into a per-session folder under their original names.
  upload_dir <- file.path(tempdir(), paste0("upload_", session$token))
  session$onSessionEnded(function() unlink(upload_dir, recursive = TRUE))
  
  observeEvent(input$dir1, {
    files <- input$dir1
    req(files, nrow(files) > 0)
    
    unlink(upload_dir, recursive = TRUE)
    dir.create(upload_dir, recursive = TRUE, showWarnings = FALSE)
    
    if (nrow(files) == 1 && grepl("\\.zip$", files$name, ignore.case = TRUE)) {
      # Normal case: the browser packed the folder into one zip
      n_ok <- length(utils::unzip(files$datapath, exdir = upload_dir, junkpaths = TRUE))
    } else {
      # Fallback: files uploaded one by one
      nm  <- basename(files$name)
      dup <- duplicated(nm)                     # same name in two sub-folders
      nm[dup] <- paste0(seq_len(sum(dup)), "_", nm[dup])
      n_ok <- sum(file.copy(files$datapath, file.path(upload_dir, nm), overwrite = TRUE))
    }
    
    dir1(upload_dir)
    output$status <- renderText(
      paste0(n_ok, " file(s) uploaded. Click 'Load data'.")
    )
  })
  
  # ---- 6.4 Load data ----
  observeEvent(input$start, {
    DT_reac(NULL)
    showNotification("Loading results...", closeButton = FALSE, duration = 2)
    
    tryCatch({
      compiled <- isTRUE(input$Compiled_F)
      
      # Path to the uploaded data (single file or folder)
      if (compiled) {
        req(input$compiled_file)
        path1 <- input$compiled_file$datapath
      } else {
        req(dir1())
        path1 <- dir1()
      }
      ###### !!! for ONLINE APP 
      # Read the tables, then fix time zones / dates / recorder names
      # While loading, R's default time zone = the recorder time zone chosen
      # by the user. Any time the loader reads without an explicit time zone
      # is then read correctly, whatever the server's own time zone is.
      old_tz <- Sys.getenv("TZ")
      Sys.setenv(TZ = input$device_tz)
      DT <- tryCatch(
        load_selection_tables_ONLINE(
          dir1 = path1, dir2 = NULL, compiled = compiled,
          device_tz = input$device_tz, deployment_tz = input$deployment_tz
        ),
        finally = Sys.setenv(TZ = old_tz)   # always put the server setting back
      )
      if (is.null(DT) || nrow(DT) == 0) {
        showNotification("No data found or tables are empty.", type = "error")
        output$status <- renderText("No data loaded.")
        return()
      }
      DT <- apply_timezones(DT, input$device_tz, input$deployment_tz, compiled)
      
      # Fill the species list
      species <- sort(unique(DT$Common.Name))
      updateSelectizeInput(session, "species_pheno",
                           choices = c("All species", species),
                           selected = "All species", server = TRUE)
      
      DT_reac(DT)
      output$status <- renderText(
        paste("Loaded:", nrow(DT), "detections |", length(species), "species")
      )
      
    }, error = function(e) {
      where <- paste(deparse(conditionCall(e)), collapse = " ")
      showNotification(paste0("Error loading data: ", e$message, "\nIn: ", where),
                       type = "error", duration = 15)
      output$status <- renderText("Error loading data.")
    })
  })
  
  # ---- 6.5 Re-apply time zones without reloading the files ----
  observeEvent(input$restart_timezone, {
    req(DT_reac(), input$device_tz, input$deployment_tz)
    tryCatch({
      DT_reac(apply_timezones(DT_reac(), input$device_tz, input$deployment_tz,
                              isTRUE(input$Compiled_F)))
      showNotification(paste("Time zones updated:", input$device_tz, "->", input$deployment_tz),
                       type = "message", duration = 3)
    }, error = function(e) {
      showNotification(paste("Error updating time zones:", e$message), type = "error")
    })
  })
  
  # ---- 6.6 Recorder selector (only shown if a 'recorder' column exists) ----
  output$recorder_ui <- renderUI({
    req(DT_reac())
    df <- DT_reac()
    if (!"recorder" %in% names(df)) return(NULL)
    recs <- sort(unique(na.omit(df$recorder)))
    if (length(recs) == 0) return(NULL)
    selectInput("selected_recorders", "Recorder(s)",
                choices = recs, selected = recs, multiple = TRUE)
  })
  
  # ---- 6.7 Phenology plot (built when "Generate plot" is clicked) ----
  pheno_data <- eventReactive(input$start_pheno, {
    req(DT_reac())
    validate(need(!is.null(input$device_tz) && !is.null(input$deployment_tz),
                  "Please select the time zones."))
    
    # Time axis: recorder clock (no DST) or local time
    plot_tz <- if (isTRUE(input$fixed_recorder_clock)) input$device_tz else input$deployment_tz
    
    dt_pheno  <- DT_reac()
    xlim_plot <- c(min(dt_pheno$Date, na.rm = TRUE) - 1,
                   max(dt_pheno$Date, na.rm = TRUE) + 1)
    
    # Filter 1: confidence
    Voc <- dt_pheno[Confidence >= input$Confid_Pheno]
    
    # Filter 2: recorders
    if ("recorder" %in% names(Voc) && length(input$selected_recorders) > 0) {
      Voc <- Voc[recorder %in% input$selected_recorders]
    }
    
    # Filter 3: species
    sp_to_plot <- input$species_pheno
    if (length(sp_to_plot) == 0) sp_to_plot <- "All species"
    if (!"All species" %in% sp_to_plot) {
      Voc <- Voc[Common.Name %in% sp_to_plot]
    }
    
    validate(need(nrow(Voc) > 0,
                  "No detections for this selection. Lower the confidence threshold or change species / recorders."))
    
    tryCatch({
      pheno_matrix(
        Voc                  = Voc,
        SP                   = sp_to_plot,
        Unit                 = input$Unit,
        sunrise              = TRUE,
        LAT                  = input$Lat,
        LONG                 = input$Lon,
        TimeZone             = plot_tz,
        xlim_plot            = xlim_plot,
        nocturnal            = input$Noctu_plot,
        tile_alpha           = input$tile_alpha,
        log_colour           = input$log_colour,
        fixed_recorder_clock = FALSE,  # DST is handled through plot_tz
        sun_as_shapes        = TRUE    # night / twilight drawn below the tiles in plotly
      )
    }, error = function(e) {
      showNotification(paste("Phenology error:", e$message), type = "error", duration = 10)
      NULL
    })
  })
  
  # Interactive plot: hover a tile to see date, time slot and number of detections
  output$pheno_plot <- renderPlotly({
    p <- pheno_data()
    
    # Nothing to plot: show the message returned by pheno_matrix (or a default one)
    validate(need(inherits(p, "gg"),
                  if (is.character(p)) p else "No data to plot. Check the parameters and click 'Generate plot'."))
    
    pl <- plotly::ggplotly(p, tooltip = "text")
    
    # Night (dark blue) and dawn / dusk (light blue) as shapes BELOW the tiles
    pl <- pheno_add_sun_shapes(pl, p)
    
    # Transparency of the detection tiles. ggplotly ignores ggplot's 'alpha' on
    # heatmaps, so it is set directly on the plotly trace. Reading the slider
    # here also updates the plot live, without clicking "Generate plot" again.
    tile_opacity <- input$tile_alpha
    
    # Tooltips: only the detection tiles; background bands and empty cells are skipped
    for (i in seq_along(pl$x$data)) {
      tr <- pl$x$data[[i]]
      if (identical(tr$type, "heatmap")) {
        pl$x$data[[i]]$opacity       <- tile_opacity
        pl$x$data[[i]]$hoverongaps   <- FALSE                  # no tooltip on empty cells
        pl$x$data[[i]]$hovertemplate <- "%{text}<extra></extra>"
      } else if (is.null(tr$text) || all(is.na(tr$text) | tr$text == "")) {
        pl$x$data[[i]]$hoverinfo <- "skip"                     # e.g. colour bar helper traces
      } else {
        pl$x$data[[i]]$opacity <- tile_opacity                 # tiles drawn as another trace type
      }
    }
    
    pl |>
      plotly::layout(
        hoverlabel  = list(bgcolor = "white", font = list(size = 13, color = "#222222")),
        annotations = list(list(
          text = "Dark blue: night  \u00b7  Light blue: dawn / dusk twilight",
          xref = "paper", yref = "paper", x = 0, y = -0.13,
          xanchor = "left", yanchor = "top", showarrow = FALSE,
          font = list(size = 12, color = "#555555")
        )),
        margin = list(b = 90)
      ) |>
      plotly::config(displaylogo = FALSE,
                     modeBarButtonsToRemove = c("lasso2d", "select2d"),
                     toImageButtonOptions = list(format = "png", filename = "phenology_plot",
                                                 width = 1400, height = 800))
  })
}


# --- 7. LAUNCH ----------------------------------------------------------------
shinyApp(ui = ui, server = server)
