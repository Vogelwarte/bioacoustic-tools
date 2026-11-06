################################################################################
# BirdNET Filtering App - ONLINE (shinyapps.io) version
#
# Purpose:
#   Interactive Shiny application to:
#   - Load BirdNET selection tables (.txt / .csv)
#   - Apply species-specific confidence thresholds
#   - Filter detections by recorder and/or date
#   - Visualize detections through interactive plots
#   - Compare recorder performance
#   - Export filtered datasets and summaries
#
# Authors:
#   Amandine Serrurier
#   Jean-Nicolas Pradervand
#   Christophe Sahli
#
# Institution:
#   Swiss Ornithological Institute
#
# Performance notes (Nov 2026 rewrite, same results as before, much faster):
#   - Recorder / date / time / GPS are parsed once per recording (unique
#     filename) instead of once per detection, then attached to the detections.
#   - Files are read only once: the "Remove NoCall" and "GPS" checkboxes no
#     longer trigger a full re-read of all files.
#   - All filtering, counting and exports use data.table (fwrite for the
#     filtered file, grouped counts for the Excel summaries).
#
################################################################################


# Remove file upload size limitation (10 GB)
options(shiny.maxRequestSize = 10 * 1024^3)

################################################################################
# 1. LOAD REQUIRED PACKAGES
#

library(shiny)
library(data.table)
library(readxl)
library(stringr)
library(plotly)
library(leaflet)
library(writexl)
library(RColorBrewer)
library(lubridate)
library(ggplot2)
library(scales)
library(bslib)
library(DT)

# Use all available cores for fread / data.table operations
data.table::setDTthreads(0)

################################################################################
# 2. DATA HELPERS (plain functions, independent of Shiny)
#

# Same column-name cleaning as before
clean_names <- function(x) tolower(gsub("[^a-zA-Z0-9]", "_", x))

# Read one or several BirdNET selection tables and stack them
read_birdnet_tables <- function(paths, names, on_file = NULL) {
  dt_list <- lapply(seq_along(paths), function(i) {
    if (!is.null(on_file)) on_file(i, names[i])
    tryCatch({
      dt <- data.table::fread(paths[i], sep = "auto", header = TRUE, showProgress = FALSE)
      data.table::setnames(dt, clean_names(names(dt)))
      dt
    }, error = function(e) {
      cat("Reading error (", names[i], "):", conditionMessage(e), "\n")
      NULL
    })
  })
  dt_list <- dt_list[!vapply(dt_list, is.null, logical(1))]
  if (length(dt_list) == 0) return(NULL)
  data.table::rbindlist(dt_list, use.names = TRUE, fill = TRUE)
}

# Parse recorder / date / time / GPS from filenames.
# Same rules as the original app (RECORDER_DATE_TIME or Lat.._Long.._RECORDER_DATE_TIME),
# but applied to each UNIQUE filename only.
parse_filename_info <- function(filenames) {
  u <- unique(filenames)
  parts <- strsplit(u, "_")
  has_gps <- vapply(parts, function(x) {
    length(x) >= 3 && grepl("^Lat", x[1]) && grepl("^Long", x[2])
  }, logical(1))
  pick <- function(i_gps, i_plain) {
    vapply(seq_along(parts), function(k) {
      x <- parts[[k]]
      i <- if (has_gps[k]) i_gps else i_plain
      if (length(x) >= i) x[i] else NA_character_
    }, character(1))
  }
  recorder <- pick(3, 1)
  date_str <- pick(4, 2)
  time_str <- gsub("\\.wav$|\\.txt$", "", pick(5, 3))

  data.table::data.table(
    filename = u,
    lat      = as.numeric(str_match(u, "Lat(-?\\d+\\.\\d+)")[, 2]),
    long     = as.numeric(str_match(u, "Long(-?\\d+\\.\\d+)")[, 2]),
    recorder = recorder,
    datetime = as.POSIXct(paste0(date_str, time_str), format = "%Y%m%d%H%M%S", tz = "UTC"),
    date     = as.Date(date_str, "%Y%m%d")
  )
}

# Add the derived columns (species names, filename, GPS, recorder, date) to the
# stacked tables. Modifies dt by reference and returns it.
prepare_birdnet <- function(dt) {
  if (is.null(dt)) return(NULL)
  if (!"common_name" %in% names(dt)) {
    cat("Error: Column 'common_name' not found. Columns:", paste(names(dt), collapse = ", "), "\n")
    return(NULL)
  }
  dt[, common_name_original := common_name]
  dt[, common_name := tolower(trimws(common_name))]
  dt[, filename := basename(begin_path)]

  info <- parse_filename_info(dt$filename)
  idx <- match(dt$filename, info$filename)
  dt[, `:=`(
    lat      = info$lat[idx],
    long     = info$long[idx],
    recorder = info$recorder[idx],
    datetime = info$datetime[idx],
    date     = info$date[idx]
  )]
  dt
}

# Species counts, sorted like table() + order(decreasing = TRUE) did before
species_counts <- function(dt) {
  cnt <- dt[!is.na(common_name_original), .N, by = common_name_original]
  res <- data.frame(Species = cnt$common_name_original, Count = cnt$N, stringsAsFactors = FALSE)
  res <- res[order(res$Species), , drop = FALSE]
  res <- res[order(res$Count, decreasing = TRUE), , drop = FALSE]
  rownames(res) <- NULL
  res
}

# Write a tab-separated file with exactly the same layout as
# write.table(df, file, sep = "\t", row.names = FALSE), but with fwrite (much faster):
# text columns and header in double quotes, dates/times unquoted, NA written as NA.
write_like_write_table <- function(df, file) {
  out <- data.table::copy(df)
  txt_cols <- names(out)[vapply(out, function(x) is.character(x) || is.factor(x), logical(1))]
  for (col in txt_cols) {
    x <- as.character(out[[col]])
    data.table::set(out, j = col,
                    value = ifelse(is.na(x), NA_character_, paste0('"', gsub('"', '\\\\"', x), '"')))
  }
  # Date-times: same text as write.table (as.character drops "00:00:00" at midnight)
  time_cols <- names(out)[vapply(out, inherits, logical(1), what = "POSIXt")]
  for (col in time_cols) data.table::set(out, j = col, value = as.character(out[[col]]))
  data.table::setnames(out, paste0('"', names(out), '"'))
  data.table::fwrite(out, file, sep = "\t", quote = FALSE, na = "NA", dateTimeAs = "write.csv")
}

# Excel workbook: total counts + one sheet per day
counts_workbook <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) {
    return(list(Total = data.frame(Species = character(), Count = integer())))
  }
  sheets <- list(Total = species_counts(dt))
  daily <- dt[!is.na(date), .(common_name_original, date)]
  for (d in split(daily, by = "date", sorted = TRUE, keep.by = TRUE)) {
    sheets[[format(d$date[1], "%Y-%m-%d")]] <- species_counts(d)
  }
  sheets
}

################################################################################
# --- 3. UI Dark theme ---
dark_theme <- bs_theme(
  version = 5,
  bootswatch = "darkly",
  bg = "#1a1a1a",
  fg = "#e0e0e0",
  primary = "#375a7f",
  secondary = "#444444"
)

ui <- fluidPage(
  theme = dark_theme,
  div(
    style = "display: flex; justify-content: space-between; align-items: center; padding: 10px 20px; background-color: #0f0f0f; border-bottom: 1px solid #333; margin-bottom: 20px;",
    div(style = "font-size: 24px; font-weight: bold; color: #fff;", "Bioacoustic tools - Check and filter output data"),
    div(tags$img(src = "logo.png", height = "70px", style = "margin-right: 10px; border-radius: 12px; padding: 6px 10px; background-color: white;"))
  ),

  sidebarLayout(
    sidebarPanel(
      fileInput("txtfiles",
                "Upload BirdNET Selection Table(s) (.txt or .csv)",
                accept = c(".txt", ".csv"),
                multiple = TRUE),
      textOutput("loaded_files_info"),
      br(),
      p("Do you have coordinates in filenames?"),
      checkboxInput("gps_mode", "Filenames contain GPS coordinates", value = FALSE),
      checkboxInput(
        "remove_nocall",
        "Remove 'NoCall' detections",
        value = TRUE
      ),
      checkboxInput(
        "hide_zero_species",
        "Hide species with 0 detections after filtering",
        value = TRUE
      ),
      h4("Step 1: Download species list for confidence threshold filtering"),
      downloadButton("downloadSpeciesTemplate", "Download xls"),

      h4("Download daily and total raw counts for all species"),
      downloadButton("downloadRawCounts", "Download raw counts"),

      h4("Step 2: Upload species list with custom confidence thresholds"),
      fileInput("xlsfile", "Upload xls/xlsx", accept = c(".xls", ".xlsx")),

      uiOutput("date_ui"),
      uiOutput("recorder_ui"),

      h4("Step 3: Download filtered data"),
      downloadButton("downloadData", "Download filtered BirdNET file"),
      downloadButton("downloadSummary", "Download filtered counts for all species")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("READ ME",
                 h3("BirdNET Filtering App - Documentation"),
                 p("This Shiny app helps you explore, filter, and visualize your BirdNET data in a simple way.") ,
                 p("With is you can:"),
                 p("📥 Upload your BirdNET selection tables.") ,
                 p("🎯 Apply custom confidence thresholds per species."),
                 p("📈 Visualize detections over time") ,
                 p("📄 Export clean, filtered results.") ,
                 p("The app works with combined selection tables from BirdNET."),
                 h4("Prerequisites"),
                 p("Ensure R packages: shiny, data.table, readxl, stringr, plotly, leaflet, writexl, RColorBrewer, lubridate, scales, bslib, ggplot2, DT"),
                 h4("Usage Instructions"),
                 tags$ol(
                   tags$li("Upload your data: Upload one or more combined BirdNET selection tables (.txt files)."),
                   tags$li("If your filenames include coordinates like Lat-XX_Long-XX_..., check the GPS box to use location data."),
                   tags$li("📥 Download the XLS template and ✍️ Add your custom confidence thresholds per species."),
                   tags$li("📤 Upload the file back into the app. 💡 Your graphs will automatically update based on these filters!"),
                   tags$li("Download XLS of raw detections for your species."),
                   tags$li("📅 Filter by date range or by 🎙️ recorder"),
                   tags$li("Download filtered results using sidebar buttons: filtered BirdNET file and filtered counts (daily and over the whole period)"),

                   tags$li("Overview: get a quick summary of your filtered data and species detections.."),
                   tags$li("Filtering impact: See how your custom thresholds affect each species."),
                   tags$li("Recorder Comparison: Compare species detections across different recorders."),
                   tags$li("Map: Display recorder locations (available when GPS mode is enabled).")
                 )
        ),
        tabPanel("Overview",
                 h4("Preview of data"),
                 DTOutput("preview"),
                 h4("Species summary"),
                 DTOutput("species_summary")
        ),
        tabPanel("Visualizations",

                 numericInput(
                   "top_species_plot",
                   "Number of most frequent species to display",
                   value = 20,
                   min = 1,
                   step = 1
                 ),

                 h4("Occurrences per species"),
                 plotlyOutput("plot_species"),

                 h4("Occurrences per date"),
                 plotlyOutput("plot_date")
        ),
        tabPanel("Filtering impact",
                 h4("Global filtering summary"),
                 verbatimTextOutput("filtering_summary"),
                 h4("Kept vs removed detections (%)"),
                 plotlyOutput("filtering_stacked_plot")
        ),
        tabPanel("Recorder Comparison",

                 numericInput(
                   "top_species_heatmap",
                   "Number of most frequent species to display",
                   value = 30,
                   min = 1,
                   step = 1
                 ),

                 h4("Species counts per recorder"),

                 plotOutput(
                   "recorder_barplot",
                   height = "600px"
                 )
        ),
        tabPanel("Map",
                 h4("Recorder positions (GPS mode only)"),
                 leafletOutput("map_positions", height = 600)
        ),
        tabPanel("Reference",
                 fluidRow(
                   column(2, textOutput("contributions")),
                   column(8, tags$p("")),
                   column(2, textOutput("license"))
                 )
        )
      )
    )
  )
)

# --- SERVER ---
server <- function(input, output, session) {

  output$contributions <- renderText("Contributions: Amandine Serrurier, Jean-Nicolas Pradervand & Christophe Sahli\nSwiss Ornithological Institute")
  output$license <- renderText("MIT License © 2026 Jean-Nicolas Pradervand")

  # Get correct colors based on number of species
  get_safe_colors <- function(n) {
    if (n <= 0) return("#cccccc")
    if (n == 1) return("#377eb8")
    if (n == 2) return(c("#377eb8", "#e41a1c"))
    if (n <= 8) return(brewer.pal(n, "Set2"))
    return(colorRampPalette(brewer.pal(8, "Set2"))(n))
  }

  # --- 1. Uploaded files ---
  selected_files <- reactiveVal(NULL)

  observeEvent(input$txtfiles, {
    selected_files(data.frame(name = input$txtfiles$name,
                              datapath = input$txtfiles$datapath,
                              stringsAsFactors = FALSE))
  })


  # --- 2. Loading BirdNET tables (read + parse ONCE per file selection) ---
  loadedBirdNET <- reactive({
    files <- selected_files()
    if (is.null(files) || nrow(files) == 0) return(NULL)

    t0 <- Sys.time()
    dt <- withProgress(message = "Reading BirdNET tables", value = 0, {
      n <- nrow(files)
      raw <- read_birdnet_tables(
        files$datapath, files$name,
        on_file = function(i, nm) incProgress(0.7 / n, detail = nm)
      )
      setProgress(0.8, detail = "Parsing recorders and dates")
      prepare_birdnet(raw)
    })
    if (!is.null(dt)) {
      cat(sprintf("Loaded %s detections from %d file(s) in %.1f s\n",
                  format(nrow(dt), big.mark = "'"), nrow(files),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
    dt
  })

  output$loaded_files_info <- renderText({
    dt <- loadedBirdNET()
    files <- selected_files()
    if (is.null(files)) return("")
    if (is.null(dt)) return("Could not read the selected file(s).")
    sprintf("%d file(s), %s detections loaded", nrow(files), format(nrow(dt), big.mark = "'"))
  })

  # Options applied on the loaded data (no re-reading of the files)
  rawBirdNET <- reactive({
    dt <- loadedBirdNET()
    if (is.null(dt)) return(NULL)

    if (isTRUE(input$remove_nocall)) {
      dt <- dt[!common_name %chin% c("nocall", "no call")]   # new table
    } else if (!isTRUE(input$gps_mode)) {
      dt <- copy(dt)                                          # keep loaded data intact
    }
    if (!isTRUE(input$gps_mode)) {
      dt[, c("lat", "long") := NA_real_]
    }
    dt
  })

  # --- 3. UI DYNAMIC ---
  # Recorder selection
  output$recorder_ui <- renderUI({
    df <- rawBirdNET()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    recs <- sort(unique(na.omit(df$recorder)))
    if(length(recs) == 0) return(NULL)
    selectInput("selected_recorders", "Select recorder(s):", choices = recs, selected = recs, multiple = TRUE)
  })

  output$date_ui <- renderUI({
    df <- rawBirdNET()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    dates <- sort(unique(df$date))
    if(length(dates) == 0) return(NULL)

    # Filter by date
    tagList(
      radioButtons("date_mode", "Filter by:", choices = c("All" = "all", "Single day" = "single", "Period" = "range"), inline = TRUE),
      conditionalPanel("input.date_mode == 'single'",
                       dateInput("selected_date", "Choose a day:", min = min(dates), max = max(dates), value = min(dates))),
      conditionalPanel("input.date_mode == 'range'",
                       dateRangeInput("selected_range", "Choose period:", start = min(dates), end = max(dates), min = min(dates), max = max(dates)))
    )
  })

  # --- 4. Filtering ---
  confidenceData <- reactive({
    if (is.null(input$xlsfile)) return(NULL)
    df_xl <- readxl::read_excel(input$xlsfile$datapath)
    names(df_xl) <- clean_names(names(df_xl))

    name_col <- grep("common|species|name", names(df_xl), ignore.case = TRUE, value = TRUE)[1]
    thresh_col <- grep("confid|threshold|seuil", names(df_xl), ignore.case = TRUE, value = TRUE)[1]

    if (is.na(name_col)) name_col <- names(df_xl)[1]
    if (is.na(thresh_col)) thresh_col <- names(df_xl)[2]

    names(df_xl)[names(df_xl) == name_col] <- "common_name"
    names(df_xl)[names(df_xl) == thresh_col] <- "Confidence_threshold"

    df_xl$common_name <- tolower(trimws(df_xl$common_name))
    return(as.data.frame(df_xl[, c("common_name", "Confidence_threshold")]))
  })

  filteredData <- reactive({
    df <- rawBirdNET()
    if (is.null(df)) return(NULL)

    keep <- rep(TRUE, nrow(df))

    df_conf <- confidenceData()
    if (!is.null(df_conf)) {
      thr <- df_conf$Confidence_threshold[match(df$common_name, df_conf$common_name)]
      keep <- keep & !is.na(thr) & df$confidence >= thr
    }

    if (!is.null(input$date_mode)) {
      if (input$date_mode == "single" && !is.null(input$selected_date)) {
        keep <- keep & df$date == input$selected_date
      } else if (input$date_mode == "range" && !is.null(input$selected_range)) {
        keep <- keep & df$date >= input$selected_range[1] & df$date <= input$selected_range[2]
      }
    }

    if (!is.null(input$selected_recorders) && length(input$selected_recorders) > 0) {
      keep <- keep & df$recorder %chin% input$selected_recorders
    }

    df[which(keep)]
  })

  # Raw counts per species (computed once, reused by several outputs)
  rawSpeciesCounts <- reactive({
    df_raw <- rawBirdNET()
    if (is.null(df_raw)) return(NULL)
    df_raw[!is.na(common_name_original), .N, by = common_name_original]
  })

  # --- 5. Summary ---
  summaryData <- reactive({
    df_raw <- rawBirdNET()
    df_filt <- filteredData()

    # Case 1 : No raw data
    if (is.null(df_raw) || nrow(df_raw) == 0) {
      return(data.frame(common_name_original = character(), Occurrence = integer()))
    }

    all_sp <- unique(df_raw$common_name_original)

    # Case 2 : Empty filtered dataset
    if (is.null(df_filt) || nrow(df_filt) == 0) {
      return(data.frame(common_name_original = all_sp, Occurrence = 0, stringsAsFactors = FALSE))
    }

    # Case 3 : counts with data.table, all species kept (even 0)
    cnt <- df_filt[!is.na(common_name_original), .N, by = common_name_original]
    final_res <- data.frame(common_name_original = all_sp, stringsAsFactors = FALSE)
    final_res$Occurrence <- cnt$N[match(final_res$common_name_original, cnt$common_name_original)]
    final_res$Occurrence[is.na(final_res$Occurrence)] <- 0L
    final_res <- final_res[order(final_res$common_name_original), , drop = FALSE]
    rownames(final_res) <- NULL

    if (input$hide_zero_species) {
      final_res <- final_res[final_res$Occurrence > 0, , drop = FALSE]
    }

    return(final_res)
  })

  # --- 6. DISPLAY TABLES ---

  output$preview <- renderDT({
    df <- rawBirdNET()

    if (is.null(df)) {
      return(datatable(data.frame(Message = "Waiting for tables..."), options = list(dom = 't')))
    }
    if (nrow(df) == 0) {
      return(datatable(data.frame(Message = "Empty file or reading error."), options = list(dom = 't')))
    }

    datatable(df, rownames = FALSE)
  })

  output$species_summary <- renderDT({
    df <- summaryData()

    if (is.null(df) || nrow(df) == 0) {
      return(datatable(data.frame(Message = "No data to summarise."), options = list(dom = 't')))
    }

    # Sort
    df <- df[order(df$Occurrence, decreasing = T), ]

    datatable(df, rownames = FALSE)
  })

  # --- 7. Graphs ---
  output$plot_species <- renderPlotly({
    df <- summaryData()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df <- df[order(-df$Occurrence), ]

    if (!is.null(input$top_species_plot)) {
      df <- head(df, input$top_species_plot)
    }
    df$common_name_original <- factor(df$common_name_original, levels = df$common_name_original)
    cols <- get_safe_colors(nrow(df))
    plot_ly(df, x = ~common_name_original, y = ~Occurrence, type = "bar", marker = list(color = cols)) %>%
      plotly::layout(xaxis = list(title = "Species", tickangle = -45, categoryorder = "trace"), yaxis = list(title = "Count"))
  })

  output$plot_date <- renderPlotly({
    df <- filteredData()
    if (is.null(df) || nrow(df) == 0) return(NULL)

    agg <- df[!is.na(date), .(Occurrence = .N), by = date][order(date)]

    plot_ly(agg, x = ~date, y = ~Occurrence, type = "bar", marker = list(color = 'darkorange')) %>%
      plotly::layout(xaxis = list(title = "Date"), yaxis = list(title = "Count"))
  })

  output$filtering_stacked_plot <- renderPlotly({
    df_raw <- rawBirdNET()
    df_filt <- filteredData()

    if (is.null(df_raw) || nrow(df_raw) == 0) return(NULL)

    raw_cnt <- rawSpeciesCounts()
    res <- data.frame(common_name_original = raw_cnt$common_name_original,
                      Raw = raw_cnt$N, stringsAsFactors = FALSE)

    if (is.null(df_filt) || nrow(df_filt) == 0) {
      res$Kept <- 0
    } else {
      filt_cnt <- df_filt[!is.na(common_name_original), .N, by = common_name_original]
      res$Kept <- filt_cnt$N[match(res$common_name_original, filt_cnt$common_name_original)]
      res$Kept[is.na(res$Kept)] <- 0
    }
    res <- res[order(res$common_name_original), , drop = FALSE]
    res$Removed <- res$Raw - res$Kept
    res$Kept_pct <- ifelse(res$Raw > 0, res$Kept/res$Raw*100, 0)
    res$Removed_pct <- ifelse(res$Raw > 0, res$Removed/res$Raw*100, 0)

    # Prepare data for Plotly (long format)
    df_long <- data.frame(
      common_name_original = rep(res$common_name_original, 2),
      Type = c(rep("Kept", nrow(res)), rep("Removed", nrow(res))),
      Percentage = c(res$Kept_pct, res$Removed_pct),
      stringsAsFactors = FALSE
    )
    df_long$Type <- factor(df_long$Type, levels = c("Kept", "Removed"))

    plot_ly(df_long, x = ~common_name_original, y = ~Percentage, color = ~Type, type = "bar") %>%
      plotly::layout(barmode = "stack", yaxis = list(range = c(0, 100)))
  })

  output$recorder_barplot <- renderPlot({
    df <- filteredData()
    if (is.null(df) || nrow(df) == 0) return(NULL)

    # 1. Aggregate (rows with missing recorder / species are ignored, as before)
    agg <- as.data.frame(
      df[!is.na(recorder) & !is.na(common_name_original) & !is.na(begin_path),
         .(Count = .N), by = .(recorder, common_name_original)]
    )
    if (nrow(agg) == 0) return(NULL)

    # 2. Total detections per species across all recorders, sorted by abundance
    total_counts <- aggregate(Count ~ common_name_original, data = agg, FUN = sum)
    total_counts <- total_counts[order(total_counts$Count, decreasing = TRUE), ]

    # Keep only Top N species
    n_top <- min(input$top_species_heatmap, nrow(total_counts))
    top_species <- total_counts$common_name_original[1:n_top]

    agg <- agg[agg$common_name_original %in% top_species, ]

    # Preserve abundance order in heatmap
    agg$common_name_original <- factor(agg$common_name_original, levels = rev(top_species))

    # 3. Graph
    ggplot(agg, aes(x = recorder, y = common_name_original, fill = Count)) +
      # Tuiles
      geom_tile(color = "#1a1a1a", alpha = 0.95) +

      # Values
      geom_text(aes(label = Count),
                color = "grey30",
                size = 5,
                fontface = "bold") +

      # Gradient
      scale_fill_gradient(low = "#f0f9e8", high = "#2ca25f", name = "Detections",
                          trans = "log10", na.value = "grey90") +

      # Duplicate  X axis
      scale_x_discrete(position = "top") +

      labs(title = "Species Detection Heatmap per Recorder",
           subtitle = "Sorted by abundance (most abundant at top)",
           x = "Recorder", y = "Species") +
      theme_minimal(base_size = 14) +
      theme(
        # --- AXis X ---

        # Configure top
        axis.title.x.top = element_text(size = 13, face = "bold", color = "#1a1a1a", margin = margin(b = 5)),
        axis.text.x.top = element_text(angle = 45, hjust = 1, vjust = 1, size = 13, color = "#1a1a1a"),
        axis.ticks.x.top = element_line(color = "grey50", linewidth = 0.5),

        # Configure bottom
        axis.title.x.bottom = element_text(size = 13, face = "bold", color = "#1a1a1a", margin = margin(t = 5)),
        axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1, size = 13, color = "#1a1a1a"),
        axis.ticks.x.bottom = element_line(color = "grey50", linewidth = 0.5),


        axis.line.x.top = element_line(color = "grey70", linewidth = 0.5),
        axis.line.x.bottom = element_line(color = "grey70", linewidth = 0.5),

        # --- AXIS Y (Species) ---
        axis.text.y = element_text(size = 11, hjust = 1, color = "#1a1a1a", face = "plain"),
        axis.title.y = element_text(size = 13, face = "bold", color = "#1a1a1a"),
        axis.line.y = element_line(color = "grey70", linewidth = 0.5),
        axis.ticks.y = element_line(color = "grey70", linewidth = 0.5),

        # --- Title and Legend ---
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "#1a1a1a"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 12),

        legend.position = "bottom",
        legend.title = element_text(size = 14, color = "#1a1a1a"),
        legend.text = element_text(size = 12, color = "#1a1a1a"),
        legend.box.margin = margin(10, 0, 0, 0),

        # --- Background and Grid ---
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "grey90", color = "#1a1a1a"),

        plot.margin = margin(10, 10, 10, 20)
      )
  }, height = function() {
    df <- filteredData()
    if (is.null(df)) return(600)
    n_sp <- min(uniqueN(df$common_name_original), input$top_species_heatmap)
    return(max(600, n_sp * 40 + 120))
  })

  # --- 8. MAP ---
  output$map_positions <- renderLeaflet({
    df <- rawBirdNET()
    if (is.null(df) || !input$gps_mode) return(NULL)
    pos <- unique(df[!is.na(lat) & !is.na(long), .(lat, long, recorder, filename)])
    if (nrow(pos) == 0) return(NULL)
    leaflet(pos) %>% addTiles() %>% addCircleMarkers(lng=~long, lat=~lat, popup=~recorder, radius=5)
  })

  # --- 9. DOWNLOADS ---
  output$downloadData <- downloadHandler(
    filename = function() "Filtered_BirdNET.txt",
    content = function(file) {
      df <- filteredData()
      if (is.null(df)) {
        write.table(data.frame(), file, sep = "\t")
      } else {
        write_like_write_table(df, file)
      }
    }
  )

  output$downloadSummary <- downloadHandler(
    filename = function() "Filtered_Counts.xlsx",
    content = function(file) {
      writexl::write_xlsx(counts_workbook(filteredData()), file)
    }
  )

  output$downloadSpeciesTemplate <- downloadHandler(
    filename = function() "Template.xlsx",
    content = function(file) {
      df <- rawBirdNET()
      if (is.null(df)) {
        writexl::write_xlsx(data.frame(common_name="", Confidence=0), file)
      } else {
        tmp <- data.frame(
          common_name = sort(unique(df$common_name)),
          Confidence = 0,
          stringsAsFactors = FALSE
        )
        writexl::write_xlsx(tmp, file)
      }
    }
  )

  output$downloadRawCounts <- downloadHandler(
    filename = function() "RawCounts.xlsx",
    content = function(file) {
      writexl::write_xlsx(counts_workbook(rawBirdNET()), file)
    }
  )
}

############################
# Run app
############################

shinyApp(ui = ui, server = server)
