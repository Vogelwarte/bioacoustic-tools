################################################################################
# BirdNET Filtering App
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
# Notes:
#   - Supports multiple BirdNET selection tables
#   - Supports recorder coordinates embedded in filenames
#   - Optimized for large BirdNET datasets
#
################################################################################

# Remove file upload size limitation (2 GB)
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
library(dplyr)
library(skimr)
library(DT)

################################################################################
# --- 2. UI Dark theme ---
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
                 p("Ensure R packages: shiny, readr, readxl, dplyr, stringr, plotly, purrr, leaflet, tidyr, writexl, janitor, RColorBrewer, lubridate, scales, bslib, ggplot2"),
                 h4("Usage Instructions"),
                 tags$ol(
                   tags$li("Upload your data: Upload one or more combined BirdNET selection tables (.txt files)"),
                   tags$li("If your filenames include coordinates like Lat-XX_Long-XX_..., check the GPS box to use location data."),
                   tags$li("📥 Download the XLS template and ✍️ Add your custom confidence thresholds per species."), 
                   tags$li("📤 Upload the file back into the app. 💡 Your graphs will automatically update based on these filters!"),
                   tags$li("Download XLS of raw detections for your species."),
                   tags$li("📅 Filter by date range or by 🎙️ recorder"),
                   tags$li("Download filtered results using sidebar buttons: filtered BirdNET file and filtered counts (daily and over the whole period)"),
                   
                   tags$li("Overview: get a quick summary of your filtered data and species detections.."),
                   tags$li("Filtering impact: See how your custom thresholds affect each species."),
                   tags$li("Recorder Comparison: Compare species detections across different recorders."),
                   #tags$li("Recorder schedule: Visualize recording time (in minutes) per recorder and across the season."),
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
  
  # --- 1. loading BirdNET dataframe ---
  rawBirdNET <- reactive({
    if (is.null(input$txtfiles)) return(NULL)
    
    cat("DEBUG: Trying to read file...\n")
    
    dt_list <- lapply(seq_along(input$txtfiles$datapath), function(i) {
      tryCatch({
        dt <- fread(
          input$txtfiles$datapath[i],
          sep = "auto",
          header = TRUE,
          data.table = FALSE
        )
        dt$source_file <- input$txtfiles$name[i]
        dt
      }, error = function(e) {
        cat("Reading error:", e$message, "\n")
        NULL
      })
    })
    
    dt_list <- dt_list[!sapply(dt_list, is.null)]
    if (length(dt_list) == 0) return(NULL)
    
    df <- do.call(rbind, dt_list)
    cat("DEBUG: File read. Number of lines:", nrow(df), "\n")
    
    # Cleaning names
    names(df) <- tolower(gsub("[^a-zA-Z0-9]", "_", names(df)))
    
    if (!"common_name" %in% names(df)) {
      cat("Error: Column 'common_name' not found Column:", paste(names(df), collapse=", "), "\n")
      return(NULL)
    }
    
    df$common_name_original <- df$common_name
    df$common_name <- tolower(trimws(df$common_name))
    if (input$remove_nocall) {
      
      df <- df[
        !tolower(trimws(df$common_name)) %in%
          c("nocall", "no call"),
      ]
      
    }
    df$filename <- basename(df$begin_path)
    
    # extract GPS data
    if (input$gps_mode) {
      df$lat <- as.numeric(str_match(df$filename, "Lat(-?\\d+\\.\\d+)")[,2])
      df$long <- as.numeric(str_match(df$filename, "Long(-?\\d+\\.\\d+)")[,2])
    } else {
      df$lat <- NA_real_
      df$long <- NA_real_
    }
    
    # Parsing (to recover date, time and recorder name based on the structure RECORDER_DATE_TIME)
    parts <- strsplit(df$filename, "_")
    df$recorder <- sapply(parts, function(x) {
      if (length(x) >= 3 && grepl("^Lat", x[1]) && grepl("^Long", x[2])) x[3] else x[1]
    })
    df$date_str <- sapply(parts, function(x) {
      if (length(x) >= 3 && grepl("^Lat", x[1]) && grepl("^Long", x[2])) x[4] else x[2]
    })
    df$time_str <- sapply(parts, function(x) {
      idx <- if (length(x) >= 3 && grepl("^Lat", x[1]) && grepl("^Long", x[2])) 5 else 3
      if (length(x) >= idx) gsub("\\.wav$|\\.txt$", "", x[idx]) else NA_character_
    })
    
    df$datetime <- as.POSIXct(paste0(df$date_str, df$time_str), format = "%Y%m%d%H%M%S", tz = "UTC")
    df$date <- as.Date(df$date_str, "%Y%m%d")
    
    # Remove useless columns
    df$source_file <- NULL
    df$date_str <- NULL
    df$time_str <- NULL
    
    cat("DEBUG: File successfully loaded.\n")
    return(df)
  })
  
  # --- 2. UI DYNAMIC ---
  # Recorder selection
  output$recorder_ui <- renderUI({
    df <- rawBirdNET()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    recs <- sort(unique(na.omit(df$recorder)))
    if(length(recs) == 0) return(NULL)
    selectInput("selected_recorders", "Select recorder(s):", choices = recs, selected = recs, multiple = TRUE)
  })
  
  #remove species with 0 data
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
  
  # --- 3. Filtering ---
  confidenceData <- reactive({
    if (is.null(input$xlsfile)) return(NULL)
    df_xl <- readxl::read_excel(input$xlsfile$datapath)
    names(df_xl) <- tolower(gsub("[^a-zA-Z0-9]", "_", names(df_xl)))
    
    name_col <- grep("common|species|name", names(df_xl), ignore.case = TRUE, value = TRUE)[1]
    thresh_col <- grep("confid|threshold|seuil", names(df_xl), ignore.case = TRUE, value = TRUE)[1]
    
    if(is.null(name_col)) name_col <- names(df_xl)[1]
    if(is.null(thresh_col)) thresh_col <- names(df_xl)[2]
    
    names(df_xl)[names(df_xl) == name_col] <- "common_name"
    names(df_xl)[names(df_xl) == thresh_col] <- "Confidence_threshold"
    
    df_xl$common_name <- tolower(trimws(df_xl$common_name))
    return(df_xl[, c("common_name", "Confidence_threshold")])
  })
  
  filteredData <- reactive({
    df <- rawBirdNET()
    if (is.null(df)) return(NULL)
    
    if (!is.null(input$xlsfile) && !is.null(confidenceData())) {
      df_conf <- confidenceData()
      df <- merge(df, df_conf, by = "common_name", all.x = TRUE)
      df <- df[!is.na(df$Confidence_threshold) & df$confidence >= df$Confidence_threshold, ]
      df$Confidence_threshold <- NULL
    }
    
    if (!is.null(input$date_mode)) {
      if (input$date_mode == "single" && !is.null(input$selected_date)) {
        df <- df[df$date == input$selected_date, ]
      } else if (input$date_mode == "range" && !is.null(input$selected_range)) {
        df <- df[df$date >= input$selected_range[1] & df$date <= input$selected_range[2], ]
      }
    }
    
    if (!is.null(input$selected_recorders) && length(input$selected_recorders) > 0) {
      df <- df[df$recorder %in% input$selected_recorders, ]
    }
    
    return(df)
  })
  
  # --- 4. Summary (Conversion explicite) ---
  summaryData <- reactive({
    df_raw <- rawBirdNET()
    df_filt <- filteredData()
    
    # Case 1 : No raw data
    if (is.null(df_raw) || nrow(df_raw) == 0) {
      return(data.frame(common_name_original = character(), Occurrence = integer()))
    }
    
    # Case 2 : Empty filtered dataset
    if (is.null(df_filt) || nrow(df_filt) == 0) {
      res <- data.frame(
        common_name_original = unique(df_raw$common_name_original),
        Occurrence = 0,
        stringsAsFactors = FALSE
      )
      return(res)
    }
    
    # Case 3 : Use table() (Base R)
    counts <- table(df_filt$common_name_original)
    
    # Convertibg table in clean dataset
    res <- data.frame(
      common_name_original = names(counts),
      Occurrence = as.integer(counts),
      stringsAsFactors = FALSE
    )
    
    # Check that all species appear in dataframe (even 0)
    all_sp <- data.frame(
      common_name_original = unique(df_raw$common_name_original), 
      stringsAsFactors = FALSE
    )
    
    final_res <- merge(all_sp, res, by = "common_name_original", all.x = TRUE)
    
    # Replace NAs (species from raw but absent in filtered) by 0
    final_res$Occurrence[is.na(final_res$Occurrence)] <- 0
    
    if (input$hide_zero_species) {
      final_res <- final_res[
        final_res$Occurrence > 0,
      ]
    }
    
    return(final_res)
  })
  # --- 5. DISPLAY TABLES ---
  
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
  
  # --- 6. Graphs ---
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
      layout(xaxis = list(title = "Species", tickangle = -45, categoryorder = "trace"), yaxis = list(title = "Count"))
  })
  
  output$plot_date <- renderPlotly({
    df <- filteredData()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    counts <- table(df$date)
    
    agg <- data.frame(
      date = as.Date(names(counts), origin = "1970-01-01"), # 
      Occurrence = as.integer(counts),
      stringsAsFactors = FALSE
    )
    
    agg <- agg[order(agg$date), ]
    
    plot_ly(agg, x = ~date, y = ~Occurrence, type = "bar", marker = list(color = 'darkorange')) %>%
      layout(xaxis = list(title = "Date"), yaxis = list(title = "Count"))
  })
  
  comparisonData <- reactive({
    df_raw <- rawBirdNET()
    df_filt <- filteredData()
    if (is.null(df_raw)) return(NULL)
    
    agg_raw <- aggregate(begin_path ~ common_name_original, data = df_raw, FUN = length)
    names(agg_raw)[2] <- "Raw"
    
    if (is.null(df_filt) || nrow(df_filt) == 0) {
      agg_filt <- data.frame(common_name_original = unique(df_raw$common_name_original), Kept = 0)
    } else {
      agg_filt <- aggregate(begin_path ~ common_name_original, data = df_filt, FUN = length)
      names(agg_filt)[2] <- "Kept"
    }
    
    res <- merge(agg_raw, agg_filt, by = "common_name_original", all.x = TRUE)
    res$Kept[is.na(res$Kept)] <- 0
    res$Removed <- res$Raw - res$Kept
    res$Kept_pct <- ifelse(res$Raw > 0, res$Kept/res$Raw*100, 0)
    res$Removed_pct <- ifelse(res$Raw > 0, res$Removed/res$Raw*100, 0)
    return(res)
  })
  
  output$filtering_stacked_plot <- renderPlotly({
    df_raw <- rawBirdNET()
    df_filt <- filteredData()
    
    if (is.null(df_raw) || nrow(df_raw) == 0) return(NULL)
    
    # Raw count
    counts_raw <- table(df_raw$common_name_original)
    agg_raw <- data.frame(
      common_name_original = names(counts_raw),
      Raw = as.integer(counts_raw),
      stringsAsFactors = FALSE
    )
    
    # Filtered count 
    if (is.null(df_filt) || nrow(df_filt) == 0) {
      agg_filt <- data.frame(
        common_name_original = unique(df_raw$common_name_original),
        Kept = 0,
        stringsAsFactors = FALSE
      )
    } else {
      counts_filt <- table(df_filt$common_name_original)
      agg_filt <- data.frame(
        common_name_original = names(counts_filt),
        Kept = as.integer(counts_filt),
        stringsAsFactors = FALSE
      )
    }
    
    # Fuse
    res <- merge(agg_raw, agg_filt, by = "common_name_original", all.x = TRUE)
    res$Kept[is.na(res$Kept)] <- 0
    res$Removed <- res$Raw - res$Kept
    res$Kept_pct <- ifelse(res$Raw > 0, res$Kept/res$Raw*100, 0)
    res$Removed_pct <- ifelse(res$Raw > 0, res$Removed/res$Raw*100, 0)
    
    # Prepare data for Plotly (Melt manuel)
    df_long <- data.frame(
      common_name_original = rep(res$common_name_original, 2),
      Type = c(rep("Kept", nrow(res)), rep("Removed", nrow(res))),
      Percentage = c(res$Kept_pct, res$Removed_pct),
      stringsAsFactors = FALSE
    )
    df_long$Type <- factor(df_long$Type, levels = c("Kept", "Removed"))
    
    plot_ly(df_long, x = ~common_name_original, y = ~Percentage, color = ~Type, type = "bar") %>%
      layout(barmode = "stack", yaxis = list(range = c(0, 100)))
  })
  
  output$recorder_barplot <- renderPlot({
    df <- filteredData()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    # 1. Aggregate
    agg <- aggregate(begin_path ~ recorder + common_name_original, data = df, FUN = length)
    names(agg)[3] <- "Count"
    if (input$hide_zero_species) {
      
      total_sp <- aggregate(
        Count ~ common_name_original,
        data = agg,
        FUN = sum
      )
      
      keep_sp <- total_sp$common_name_original[
        total_sp$Count > 0
      ]
      
      agg <- agg[
        agg$common_name_original %in% keep_sp,
      ]
    }
    
    # 2. Sort by abundance decreasing
    # Total detections per species across all recorders
    total_counts <- aggregate(
      Count ~ common_name_original,
      data = agg,
      FUN = sum
    )
    
    # Sort by abundance
    total_counts <- total_counts[
      order(total_counts$Count, decreasing = TRUE),
    ]
    
    # Keep only Top N species
    n_top <- min(
      input$top_species_heatmap,
      nrow(total_counts)
    )
    
    top_species <- total_counts$common_name_original[1:n_top]
    
    # Filter heatmap dataset
    agg <- agg[
      agg$common_name_original %in% top_species,
    ]
    
    # Preserve abundance order in heatmap
    agg$common_name_original <- factor(
      agg$common_name_original,
      levels = rev(top_species)
    )
    
    # 3. Calculate dynamic height
    n_species <- length(top_species)
    plot_height <- max(600, n_species * 40) 
    
    # 4. Graph
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
        axis.ticks.x.top = element_line(color = "grey50", size = 0.5), # Ligne et ticks visibles
        
        # Configure bottom
        axis.title.x.bottom = element_text(size = 13, face = "bold", color = "#1a1a1a", margin = margin(t = 5)),
        axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1, size = 13, color = "#1a1a1a"),
        axis.ticks.x.bottom = element_line(color = "grey50", size = 0.5),
        
        
        axis.line.x.top = element_line(color = "grey70", size = 0.5),
        axis.line.x.bottom = element_line(color = "grey70", size = 0.5),
        
        # --- AXIS Y (Species) ---
        axis.text.y = element_text(size = 11, hjust = 1, color = "#1a1a1a", face = "plain"),
        axis.title.y = element_text(size = 13, face = "bold", color = "#1a1a1a"),
        axis.line.y = element_line(color = "grey70", size = 0.5),
        axis.ticks.y = element_line(color = "grey70", size = 0.5),
        
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
    n_sp <- min(
      length(unique(df$common_name_original)),
      input$top_species_heatmap
    )
    return(max(600, n_sp * 40 + 120)) 
  })
  
  
  # --- 7. SCHEDULE & MAP (Simplified) ---
  scheduleData <- reactive({
    if (is.null(input$schedule_file)) return(NULL)
    lines <- readLines(input$schedule_file$datapath)
    df <- data.frame(path = lines, stringsAsFactors = FALSE)
    df$filename <- basename(df$path)
    parts <- strsplit(df$filename, "_")
    df$recorder <- sapply(parts, `[`, 1)
    df$date_str <- sapply(parts, `[`, 2)
    df$time_str <- sapply(parts, function(x) gsub("\\.wav$|\\.WAV$", "", x[3]))
    df$datetime <- as.POSIXct(paste0(df$date_str, df$time_str), format = "%Y%m%d%H%M%S", tz = "UTC")
    df$date <- as.Date(df$date_str, "%Y%m%d")
    df[!is.na(df$date), ]
  })
  
  output$recorder_daily_plot <- renderPlot({
    df <- scheduleData()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$hour <- lubridate::hour(df$datetime)
    df$half_day <- ifelse(df$hour < 12, "00–12h", "12–24h")
    agg <- aggregate(n ~ recorder + date + half_day, data = df, FUN = length)
    agg$recording <- ifelse(agg$n > 0, "yes", "no")
    
    # Expansion manuelle
    all_comb <- expand.grid(recorder = unique(agg$recorder), date = unique(agg$date), half_day = c("00–12h", "12–24h"))
    df_plot <- merge(all_comb, agg, by = c("recorder", "date", "half_day"), all.x = TRUE)
    df_plot$recording[is.na(df_plot$recording)] <- "no"
    
    ggplot(df_plot, aes(x = date, y = half_day, fill = recording)) +
      geom_tile(color = "black") + scale_fill_manual(values = c("yes" = "#2E7D32", "no" = "#CFEFCF")) +
      facet_wrap(~ recorder, ncol = 1, strip.position = "left") + theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$schedule_plot_total <- renderPlot({
    df <- scheduleData()
    if (is.null(df) || nrow(df) == 0 || is.null(input$rec_length)) return(NULL)
    df$hour <- lubridate::hour(df$datetime)
    df$half_day <- ifelse(df$hour < 12, "00–12h", "12–24h")
    agg <- aggregate(minutes_recorded ~ date + half_day, data = df, FUN = function(x) length(x) * input$rec_length)
    names(agg)[3] <- "minutes_recorded"
    
    ggplot(agg, aes(x = date, y = half_day, fill = minutes_recorded)) +
      geom_tile() + scale_fill_gradient(low = "white", high = "darkred") + theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$map_positions <- renderLeaflet({
    df <- rawBirdNET()
    if (is.null(df) || !input$gps_mode) return(NULL)
    pos <- df[!is.na(df$lat) & !is.na(df$long), c("lat", "long", "recorder", "filename")]
    pos <- unique(pos)
    if (nrow(pos) == 0) return(NULL)
    leaflet(pos) %>% addTiles() %>% addCircleMarkers(lng=~long, lat=~lat, popup=~recorder, radius=5)
  })
  
  # --- 8. DOWNLOADS (Simplified) ---
  output$downloadData <- downloadHandler(
    filename = function() "Filtered_BirdNET.txt",
    content = function(file) {
      df <- filteredData()
      if (is.null(df)) write.table(data.frame(), file, sep="\t")
      else write.table(df, file, sep="\t", row.names=FALSE)
    }
  )
  
  output$downloadSummary <- downloadHandler(
    filename = function() "Summary.xlsx",
    content = function(file) {
      df <- summaryData()
      writexl::write_xlsx(df, file)
    }
  )
  
  output$downloadSpeciesTemplate <- downloadHandler(
    filename = function() "Template.xlsx",
    content = function(file) {
      df <- rawBirdNET()
      if (is.null(df)) {
        writexl::write_xlsx(data.frame(common_name="", Confidence=0), file)
      } else {
        tmp <- unique(df[, c("common_name")])
        tmp$Confidence <- 0
        writexl::write_xlsx(tmp, file)
      }
    }
  )
  
  output$downloadRawCounts <- downloadHandler(
    filename = function() "RawCounts.xlsx",
    content = function(file) {
      df <- rawBirdNET()
      if (is.null(df)) {
        writexl::write_xlsx(data.frame(), file)
      } else {
        agg <- aggregate(Count ~ common_name_original, data = df, FUN = length)
        names(agg)[2] <- "Total"
        writexl::write_xlsx(list(Total = agg), file)
      }
    }
  )
}

############################
# Run app
############################

shinyApp(ui = ui, server = server)