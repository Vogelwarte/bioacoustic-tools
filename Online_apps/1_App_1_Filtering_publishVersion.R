# ==============================================================================



# Remove file upload size limitation (10 GB)
options(shiny.maxRequestSize = 10 * 1024^3)

library(shiny)
library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(plotly)
library(purrr)
library(leaflet)
library(tidyr)
library(writexl)
library(janitor)
library(RColorBrewer)
library(lubridate)
library(ggplot2)
library(scales)
library(bslib)



############################
# UI
############################
# Définition du thème sombre avec bslib
dark_theme <- bs_theme(
  version = 5,
  bootswatch = "darkly", # Thème sombre prédéfini
  bg = "#1a1a1a",        # Fond principal
  fg = "#e0e0e0",        # Couleur du texte
  primary = "#375a7f",   # Couleur primaire (bleu sombre)
  secondary = "#444444"  # Couleur secondaire
)

ui <- fluidPage(
  theme = dark_theme, 
   #titlePanel("BirdNET Filtering App"),
  # --- En-tête personnalisé avec Logo ---
  div(
    style = "display: flex; justify-content: space-between; align-items: center; padding: 10px 20px; background-color: #0f0f0f; border-bottom: 1px solid #333; margin-bottom: 20px;",
    div(
      style = "font-size: 24px; font-weight: bold; color: #fff;",
      "BirdNET-ResChecker - Check and filter output data"
    ),
    div(
      tags$img(
        src = "logo.png", 
        height = "70px", 
        # padding crée l'espace, background-color remplit cet espace
        style = "margin-right: 10px; border-radius: 12px; padding: 6px 10px; background-color: white;"
      ) 
    )
  ),
  

  
  sidebarLayout(
    sidebarPanel(
      #p("Upload BirdNET combined selection table(s)."),
      fileInput("txtfiles", "Step 1: Upload BirdNET Selection Table(s) (.txt)", accept = ".txt", multiple = TRUE),
      
      p("Do you have coordinates in filenames?"),
      checkboxInput("gps_mode", "Filenames contain GPS coordinates", value = FALSE),
      
      h4("Step 2.1: Download species list for confidence threshold filtering"),
      downloadButton("downloadSpeciesTemplate", "Download xls"),
      
      h4("Step 2.2: Download daily and total raw counts for all species"),
      downloadButton("downloadRawCounts", "Download raw counts"),
      
      h4("Step 3: Upload species list with custom confidence thresholds"),
      fileInput("xlsfile", "Upload xls/xlsx", accept = c(".xls", ".xlsx")),
      
      uiOutput("date_ui"),
      uiOutput("recorder_ui"),
      
      h4("Step 4: Download filtered data"),
      downloadButton("downloadData", "Download filtered BirdNET file"),
      downloadButton("downloadSummary", "Download filtered counts for all species")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("READ ME",
                 h3("BirdNET Filtering App - Documentation"),
                 p("This Shiny app helps you explore, filter, and visualize your BirdNET data in a simple way.") ,
                 p("With this App you can:"),
                 p("📥 Upload your BirdNET selection tables.") ,
                 p("🎯 Apply custom confidence thresholds per species."),
                 p("📈 Visualize detections over time") ,
                 p("📄 Export clean, filtered results.") ,
                 p("The app works with combined selection tables from BirdNET."),
                 h4("Prerequisites"),
                 p("Ensure R packages: shiny, readr, readxl, dplyr, stringr, plotly, purrr, leaflet, tidyr, writexl, janitor, RColorBrewer, lubridate, scales, bslib, ggplot2"),
                 h4("Usage Instructions"),
                 tags$ol(
                   tags$li("Step 1.: Upload your data: Upload one or more combined BirdNET selection tables (.txt files). For this use the Browse button on the top"),
                   tags$li("If your filenames include coordinates like Lat-XX_Long-XX_..., check the GPS box to use location data.You should be able to visalize your data on a map."),
                   tags$li("📥 Step 2.1: Download an XLS file providing a list of all your species. ✍️ You can add your custom confidence thresholds per species in this template and use it to filter your dataset.For unverified data, either use 0 (all) or 1.1 (none)."), 
                   tags$li("📥 Step 2.2: Download an XLS file with a list of detections per species."), 
                   tags$li("📤 Step 3. :Upload the file back into the app. 💡 Your graphs will automatically update based on these filters!"),
                   tags$li("Step 4: Download XLS of raw detections for your species."),
                   tags$li("Options: 📅 Filter by date range or by 🎙️ recorder"),
                   tags$li("Download filtered results using sidebar buttons: filtered BirdNET file and filtered counts (daily and over the whole period)"),
                   
                   tags$li("Overview: get a quick summary of your filtered data and species detections.."),
                   tags$li("Filtering impact: See how your custom thresholds affect each species."),
                   tags$li("Recorder Comparison: Compare species detections across different recorders."),
                   tags$li("Recorder schedule: Visualize recording time (in minutes) per recorder and across the season."),
                   tags$li("Map: Display recorder locations (available when GPS mode is enabled).")
                 )
        ),
        tabPanel("Overview",
                 h4("Preview of data"),
                 tableOutput("preview"),
                 h4("Species summary"),
                 tableOutput("species_summary")
        ),
        tabPanel("Visualizations",
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
                 h4("Species counts per recorder"),
                 plotlyOutput("recorder_barplot", height = "600px")
        ),
        tabPanel("Recording Schedule",
                 numericInput("rec_length", 
                              "Define recording length (in minutes):", 
                              value = 1, min = 1, step = 1),
                 
                 h4("Upload recording file list (.txt)"),
                 fileInput("schedule_file", "Upload TXT file", accept = ".txt"),
                 
                 h4("Daily recording (hours per recorder)"),
                 plotOutput("recorder_daily_plot", height = "400px"),
                 
                 h4("Total recording density (all recorders combined)"),
                 plotOutput("schedule_plot_total", height = "400px")
        ),
        tabPanel("Map",
                 h4("Recorder positions (GPS mode only)"),
                 leafletOutput("map_positions", height = 600)
        ), 
        nav_panel("Reference",
                  layout_columns(
                    col_widths = c(2, 8, 2),
                    textOutput("contributions"),
                    textOutput("license")
                  )
        ),
      )
    )
  )
)

############################
# SERVER
############################
server <- function(input, output, session) {
  
  # --- Textes statiques ---
  output$contributions <- renderText("Contributions: Amandine Serrurier, Jean-Nicolas Pradervand & Christophe Sahli\nSwiss Ornithological Institute")
  output$license <- renderText("MIT License © 2026 Jean-Nicolas Pradervand")
  
  
  ############################
  # Load raw BirdNET data
  ############################
  rawBirdNET <- reactive({
    req(input$txtfiles)
    
    map2_dfr(input$txtfiles$datapath, input$txtfiles$name, 
             ~ read_delim(.x, delim = "\t", col_types = cols()) %>%
               mutate(source_file = .y)) %>%
      janitor::clean_names() %>%
      mutate(
        common_name_original = common_name,
        common_name = str_trim(str_to_lower(common_name)),
        
        filename = basename(begin_path),
        
        lat = parse_number(str_extract(filename, "Lat-?\\d+\\.\\d+")),
        long = parse_number(str_extract(filename, "Long-?\\d+\\.\\d+")),
        
        # parsing
        parts = str_split(filename, "_"),
        
        parts_clean = map(parts, function(x) {
          if (length(x) >= 3 && str_detect(x[1], "^Lat") && str_detect(x[2], "^Long")) {
            x[-c(1,2)]
          } else {
            x
          }
        }),
        
        recorder = as.factor(map_chr(parts_clean, ~ if(length(.x) >= 1) .x[1] else NA_character_)),
        date_str = map_chr(parts_clean, ~ if(length(.x) >= 2) .x[2] else NA_character_),
        time_str = map_chr(parts_clean, ~ if(length(.x) >= 3) str_remove(.x[3], "\\.wav$|\\.txt$") else NA_character_),
        
        datetime = as.POSIXct(
          paste0(date_str, time_str),
          format = "%Y%m%d%H%M%S",
          tz = "UTC"
        ),
        
        date = as.Date(date_str, "%Y%m%d")
      )
  }) 
  
  ############################
  # Recorder selection
  ############################
  output$recorder_ui <- renderUI({
    req(rawBirdNET())
    recs <- sort(unique(na.omit(rawBirdNET()$recorder)))
    selectInput("selected_recorders", "Select recorder(s):", choices = recs, selected = recs, multiple = TRUE)
  })
  
  ############################
  # Date selection
  ############################
  output$date_ui <- renderUI({
    req(rawBirdNET())              # only proceed if data is loaded
    dates <- sort(unique(rawBirdNET()$date))
    req(length(dates) > 0)         # only proceed if there are dates
    
    tagList(
      radioButtons("date_mode", "Filter by:", 
                   choices = c("All" = "all", "Single day" = "single", "Period" = "range"), 
                   inline = TRUE),
      conditionalPanel(
        condition = "input.date_mode == 'single'",
        dateInput("selected_date", "Choose a day:", 
                  min = min(dates), max = max(dates), value = min(dates))
      ),
      conditionalPanel(
        condition = "input.date_mode == 'range'",
        dateRangeInput("selected_range", "Choose period:", 
                       start = min(dates), end = max(dates),
                       min = min(dates), max = max(dates))
      )
    )
  })
  
  ############################
  # Confidence thresholds
  ############################
  confidenceData <- reactive({
    req(input$xlsfile)
    read_excel(input$xlsfile$datapath) %>%
      rename(Confidence_threshold = Confidence) %>%
      select(common_name, Confidence_threshold) %>%
      mutate(common_name = str_trim(str_to_lower(common_name)))
  })
  
  ############################
  # Filtered data reactive
  ############################
  filteredData <- reactive({
    req(rawBirdNET())
    df <- rawBirdNET()
    
    # Apply confidence thresholds if xls is uploaded
    if (!is.null(input$xlsfile)) {
      df <- df %>%
        left_join(confidenceData(), by = "common_name") %>%
        filter(!is.na(Confidence_threshold) & confidence >= Confidence_threshold)
    }
    
    # Apply date filtering
    if (!is.null(input$date_mode)) {
      if (input$date_mode == "single" && !is.null(input$selected_date)) {
        df <- df %>% filter(date == input$selected_date)
      } else if (input$date_mode == "range" && !is.null(input$selected_range)) {
        df <- df %>% filter(date >= input$selected_range[1], date <= input$selected_range[2])
      }
    }
    
    # Apply recorder filtering
    if (!is.null(input$selected_recorders) && length(input$selected_recorders) > 0) {
      df <- df %>% filter(recorder %in% input$selected_recorders)
    }
    
    df
  })
  
  ############################
  # Overview tab: preview & species summary
  ############################
  output$preview <- renderTable({
    df <- rawBirdNET()
    req(nrow(df) > 0)
    
    df %>%
      mutate(across(everything(), as.character)) %>%  # force safe display
      head()
  })
  
  output$species_summary <- renderTable({
    rawBirdNET() %>% 
      group_by(common_name_original) %>% 
      summarise(Occurrence = n(), .groups = "drop")
  })
  
  summaryData <- reactive({
    df <- filteredData()
    df_summary <- df %>%
      group_by(common_name_original) %>%
      summarise(Occurrence = n(), .groups = "drop")
    
    # Include all species even if removed after filtering
    all_species <- tibble(common_name_original = unique(rawBirdNET()$common_name_original))
    all_species %>%
      left_join(df_summary, by = "common_name_original") %>%
      mutate(Occurrence = ifelse(is.na(Occurrence), 0, Occurrence))
  })
  
  ############################
  # Visualizations: species & date
  ############################
  output$plot_species <- renderPlotly({
    df <- summaryData()
    if (nrow(df) == 0) return(NULL)
    
    species_list <- unique(df$common_name_original)
    base_colors <- RColorBrewer::brewer.pal(8, "Set2")
    colors <- colorRampPalette(base_colors)(length(species_list))
    
    plot_ly(df, x = ~common_name_original, y = ~Occurrence, type = "bar",
            marker = list(color = rep(colors, length.out = nrow(df)))) %>%
      layout(title = "Occurrences per species",
             xaxis = list(title = "Species", tickangle = -45),
             yaxis = list(title = "Count"))
  })
  
  output$plot_date <- renderPlotly({
    df <- filteredData()
    if (nrow(df) == 0) return(NULL)
    
    df_summary <- df %>% group_by(date) %>% summarise(Occurrence = n(), .groups = "drop")
    plot_ly(df_summary, x = ~date, y = ~Occurrence, type = "bar", marker = list(color = 'darkorange')) %>%
      layout(title = "Occurrences per date", xaxis = list(title = "Date"), yaxis = list(title = "Count"))
  })
  
  
  ############################
  # Filtering impact
  ############################
  comparisonData <- reactive({
    raw <- rawBirdNET()
    filt <- filteredData()
    
    raw_counts <- raw %>% group_by(common_name_original) %>% summarise(Raw = n(), .groups="drop")
    filt_counts <- filt %>% group_by(common_name_original) %>% summarise(Kept = n(), .groups="drop")
    
    raw_counts %>% left_join(filt_counts, by="common_name_original") %>%
      mutate(Kept = ifelse(is.na(Kept), 0, Kept),
             Removed = Raw - Kept,
             Kept_pct = ifelse(Raw>0, Kept/Raw*100, 0),
             Removed_pct = ifelse(Raw>0, Removed/Raw*100, 0))
  })
  
  output$filtering_stacked_plot <- renderPlotly({
    df <- comparisonData()
    df_long <- df %>% select(common_name_original, Kept_pct, Removed_pct) %>%
      pivot_longer(cols = c(Kept_pct, Removed_pct), names_to = "Type", values_to = "Percentage")
    
    plot_ly(df_long, x = ~common_name_original, y = ~Percentage, color = ~Type, type="bar") %>%
      layout(barmode = "stack", yaxis = list(title="Percentage", range=c(0,100)),
             xaxis = list(title="Species", tickangle=-45), title="Filtering impact (%)")
  })
  
  output$filtering_summary <- renderText({
    raw_total <- nrow(rawBirdNET())
    kept_total <- nrow(filteredData())
    removed_total <- raw_total - kept_total
    paste0("Total detections: ", raw_total, "\n",
           "Kept: ", kept_total, " (", round(kept_total/raw_total*100,1), "%)\n",
           "Removed: ", removed_total, " (", round(removed_total/raw_total*100,1), "%)")
  })
  
  ###########################
  # Species per recorders
  ###########################
  
  ## Recorder Comparison: counts per recorder without filtering
  
  
  output$recorder_barplot <- renderPlotly({
    df <- rawBirdNET() %>%
      mutate(recorder = as.factor(recorder)) %>%   
      group_by(recorder, common_name_original) %>%
      summarise(Count = n(), .groups = "drop")

    
    if(nrow(df) == 0) return(NULL)
    
    n_colors <- min(length(unique(df$common_name_original)), 8)
    colors <- RColorBrewer::brewer.pal(n_colors, "Set2")
    
    plot_ly(df, x = ~recorder, y = ~Count, color = ~common_name_original, type = "bar",
            marker = list(color = rep(colors, length.out = nrow(df)))) %>%
      layout(barmode = "stack",
             title = "Species counts per recorder",
             xaxis = list(title = "Recorder"),
             yaxis = list(title = "Count"))
  })
  
  ############################
  # Recorder comparison
  ############################
  # Single recorder data
  output$recorder_daily_plot <- renderPlot({
    df <- scheduleData()
    req(nrow(df) > 0)
    
    # Extract hour and define the two half-day bins
    df <- df %>%
      mutate(
        hour = lubridate::hour(datetime),
        half_day = ifelse(hour < 12, "00–12h", "12–24h")
      )
    
    # Compute presence/absence per recorder x date x half-day
    df_summary <- df %>%
      group_by(recorder, date, half_day) %>%
      summarise(n = n(), .groups = "drop") %>%
      mutate(recording = ifelse(n > 0, "yes", "no"))
    
    # Make sure both bins exist even when empty
    all_combinations <- expand.grid(
      recorder = unique(df_summary$recorder),
      date     = unique(df_summary$date),
      half_day = c("00–12h", "12–24h"),
      stringsAsFactors = FALSE
    )
    
    df_plot <- all_combinations %>%
      left_join(df_summary, by = c("recorder", "date", "half_day")) %>%
      mutate(recording = ifelse(is.na(recording), "no", recording))
    
    ggplot(df_plot, aes(x = date, y = half_day, fill = recording)) +
      geom_tile(color = "black") +
      scale_fill_manual(values = c("yes" = "#2E7D32", "no" = "#CFEFCF"), guide = "none") +
      facet_wrap(~ recorder, ncol = 1, strip.position = "left") +
      scale_y_discrete(limits = c("12–24h", "00–12h")) +
      labs(
        title = "Daily recording schedule per recorder",
        x = "Date",
        y = NULL
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid = element_blank(),
        axis.text.y = element_blank(),       # remove tick labels
        axis.ticks.y = element_blank(),      # remove ticks
        strip.text.y.left = element_text(size = 8, face = "bold", angle = 0),  # recorder names --> Angle is still not good!
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  })
  
# Summary for all recorders  
  
  output$schedule_plot_total <- renderPlot({
    df <- scheduleData()
    req(nrow(df) > 0)
    req(input$rec_length)
    
    rec_length <- input$rec_length
    
    df <- df %>%
      mutate(
        hour = lubridate::hour(datetime),
        half_day = ifelse(hour < 12, "00–12h", "12–24h")
      )
    
    df_summary <- df %>%
      group_by(date, half_day) %>%
      summarise(
        minutes_recorded = n() * rec_length,
        .groups = "drop"
      )
    
    ggplot(df_summary, aes(x = date, y = half_day, fill = minutes_recorded)) +
      geom_tile(color = "grey80") +
      scale_fill_gradient(
        low = "white",
        high = "darkred",
        name = "Minutes"
      ) +
      labs(
        title = "Total recording density (all recorders, 12h bins)",
        x = "Date",
        y = "Time block"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  }) 
  
    
  ############################
  # Map
  ############################
  output$map_positions <- renderLeaflet({
    req(input$gps_mode, rawBirdNET())
    positions <- rawBirdNET() %>% filter(!is.na(lat) & !is.na(long)) %>% distinct(lat, long, recorder, filename)
    if(nrow(positions)==0) return(NULL)
    
    leaflet(positions) %>% addTiles() %>%
      addCircleMarkers(lng=~long, lat=~lat, label=~as.character(recorder),
                       popup = ~paste0("<b>Recorder:</b> ", recorder, "<br>",
                                       "<b>Lat:</b> ", lat, "<br>",
                                       "<b>Long:</b> ", long, "<br>",
                                       "<b>File:</b> ", filename),
                       radius=6, color="blue", fillOpacity=0.7) %>%
      fitBounds(lng1=min(positions$long), lat1=min(positions$lat),
                lng2=max(positions$long), lat2=max(positions$lat))
  })
  ############################
  # Build schedule data
  ############################
  
  scheduleData <- reactive({
    req(input$schedule_file)
    
    readr::read_lines(input$schedule_file$datapath) %>%
      tibble(path = .) %>%
      mutate(
        filename = basename(path),
        
        # Split by underscore
        parts = str_split(filename, "_"),
        
        # Recorder: always part 1
        recorder = map_chr(parts, ~ .x[1]),
        recorder = as.factor(recorder),
        
        # Date: always part 2
        date_str = map_chr(parts, ~ .x[2]),
        
        # Time: part 3, remove extension
        time_str = map_chr(parts, ~ str_remove(.x[3], "\\.wav$|\\.WAV$")),
        
        # Build datetime
        datetime = as.POSIXct(
          paste0(date_str, time_str),
          format = "%Y%m%d%H%M%S",
          tz = "UTC"
        ),
        
        date = as.Date(date_str, "%Y%m%d")
      ) %>%
      filter(!is.na(date), !is.na(datetime))
  })
  
  
 
  
  
  ############################
  # Downloads
  ############################
  output$downloadRawCounts <- downloadHandler(
    filename = function() {
      paste0("BirdNET_raw_counts_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      df <- rawBirdNET()
      all_species <- tibble(common_name_original = unique(df$common_name_original))
      
      # Total counts per species
      all_days <- df %>%
        group_by(common_name_original) %>%
        summarise(Count = n(), .groups = "drop") %>%
        right_join(all_species, by = "common_name_original") %>%
        mutate(Count = ifelse(is.na(Count), 0, Count))
      
      # Counts per day
      days <- sort(unique(df$date))
      per_day_list <- map(days, function(d) {
        df_day <- df %>%
          filter(date == d) %>%
          group_by(common_name_original) %>%
          summarise(Count = n(), .groups = "drop") %>%
          right_join(all_species, by = "common_name_original") %>%
          mutate(Count = ifelse(is.na(Count), 0, Count))
        df_day
      })
      names(per_day_list) <- as.character(days)
      
      # Combine into sheets
      sheets <- c(list(All_Days = all_days), per_day_list)
      writexl::write_xlsx(sheets, path = file)
    }
  )
  
  
  
  output$downloadSpeciesTemplate <- downloadHandler(
    filename = function(){ paste0("BirdNET_species_template_", Sys.Date(), ".xlsx") },
    content = function(file){
      df <- rawBirdNET() %>% group_by(common_name) %>% summarise(Detections=n(), .groups="drop") %>%
        select(common_name) %>% mutate(Confidence=0)
      writexl::write_xlsx(df, file)
    }
  )
  
  output$downloadData <- downloadHandler(
    filename = function(){ "Filtered_BirdNET_selection_table.txt" },
    content = function(file){ write_delim(filteredData(), file, delim="\t") }
  )
  
  output$downloadSummary <- downloadHandler(
    filename = function(){
      recs <- paste(input$selected_recorders, collapse="_")
      paste0("Filtered_Occurrences_", recs, "_", Sys.Date(), ".xlsx")
    },
    content = function(file){
      df <- filteredData()
      all_species <- tibble(common_name_original=unique(rawBirdNET()$common_name_original))
      
      all_days <- df %>% group_by(common_name_original) %>% summarise(Count=n(), .groups="drop") %>%
        right_join(all_species, by="common_name_original") %>% mutate(Count=ifelse(is.na(Count),0,Count))
      
      days <- sort(unique(df$date))
      per_day_list <- map(days, function(d){
        df_day <- df %>% filter(date==d) %>% group_by(common_name_original) %>%
          summarise(Count=n(), .groups="drop") %>% right_join(all_species, by="common_name_original") %>%
          mutate(Count=ifelse(is.na(Count),0,Count))
        df_day
      })
      names(per_day_list) <- as.character(days)
      sheets <- c(list(All_Days=all_days), per_day_list)
      writexl::write_xlsx(sheets, path=file)
    }
  )
  
}

############################
# Run app
############################
shinyApp(ui = ui, server = server)

## 
