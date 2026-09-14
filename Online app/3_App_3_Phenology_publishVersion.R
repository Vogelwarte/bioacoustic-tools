# ==============================================================================
# BirdNET-ResChecker - PHENOLOGY
# Version Optimisée & Corrigée (Thème Sombre + Logo)
# ==============================================================================
options(shiny.maxRequestSize = 10 * 1024^3)
rsconnect::setAccountInfo(name='vogelwarte',
                          token='5FA01481ECE5F33483C37F69731D3CF9',
                          secret='bBVKSYJlP/1Xz2Ynz8FuwxotsbgoQQTRrKekRXh7')



# --- 1. CHARGEMENT DES LIBRAIRIES ---------------------------------------------
# library(pacman)
library(shiny)
library(grid)
library(bslib)
library(stringr)
library(data.table)
library(shinyFiles)
library(tidyr)
library(stringi)
library(fs)
library(parallel)
library(ggplot2)
library(mgcv)
library(gratia)
library(hms)
library(suncalc)
library(lubridate)
library(DT)
library(seewave)
library(tuneR)
library(base64enc)
library(av)
library(audio)
library(dipsaus)
library(future)

# --- 2. CHARGEMENT DES FONCTIONS PERSONNALISÉES -------------------------------
  source("Common_functions/get_roots.R")
  source("Common_functions/pheno_matrix.r")
  source("Common_functions/clean_text.R")
  source("Common_functions/export_selected_audio.R")
  source("Common_functions/prepare_for_spectro.R")
  source("Common_functions/load_selection_tables_ONLINE.R")
#test
# --- 3. CONFIGURATION GLOBALE -------------------------------------------------
# roots_home <- get_roots()
utc_timezones <- OlsonNames()

Files_message_1 <- "-- There are no results files (~.BirdNET.selection.table.txt) in the results folder. --"

# Définition du thème sombre avec bslib
dark_theme <- bs_theme(
  version = 5,
  bootswatch = "darkly", # Thème sombre prédéfini
  bg = "#1a1a1a",        # Fond principal
  fg = "#e0e0e0",        # Couleur du texte
  primary = "#375a7f",   # Couleur primaire (bleu sombre)
  secondary = "#444444"  # Couleur secondaire
)

ui <- page_sidebar( # Remplace page_fluid par page_sidebar
  theme = dark_theme,
  sidebar = sidebar(
    title = "Data Source",
    checkboxInput(
      "Compiled_F",
      "Results in compiled format?",
      value = FALSE
    ),
    
    conditionalPanel(
      condition = "input.Compiled_F == false",
      
      dipsaus::fancyDirectoryInput(
        "dir1",
        "Choose folder",
        autoCleanup = TRUE,
        autoCleanupLocked = TRUE
      )
    ),
    
    conditionalPanel(
      condition = "input.Compiled_F == true",
      
      fileInput(
        "compiled_file",
        "Choose compiled result file",
        multiple = FALSE,
        accept = c(".csv", ".txt")
      )
    ),
    
    
    # checkboxInput(
    #   "compiled_format",
    #   "Results in compiled format?",
    #   value = FALSE
    # ),
    # 
    # conditionalPanel(
    #   condition = "input.compiled_format == false",
    #   
    #   dipsaus::fancyDirectoryInput(
    #     "dir1",
    #     "Choose folder",
    #     autoCleanup = TRUE,
    #     autoCleanupLocked = TRUE
    #   )
    # ),
    # 
    # conditionalPanel(
    #   condition = "input.compiled_format == true",
    #   
    #   fileInput(
    #    "Compiled_F",
    #     "Choose compiled result file",
    #     multiple = FALSE,
    #     accept = c(
    #       ".csv",
    #       ".txt"
    #     )
    #   )
    # ),
    # dipsaus::fancyDirectoryInput(
    #   "dir1",
    #   "Select BirdNET Results Folder",
    #   autoCleanup = TRUE,
    #   autoCleanupLocked = TRUE
    # ),
    
    # shinyDirButton("dir1", "Select BirdNET Results Folder", title = "Choose Results Folder"),
    # checkboxInput("Compiled_F", label = "Results in compiled format?", value = FALSE),
    textOutput("dir1_path"),
    hr(),
    actionButton("start", "Start App", class = "btn-primary", width = "100%"),
    hr(),
    selectInput(inputId = "UTC_choice", label = "Timezone (UTC)", 
                choices = utc_timezones, selected = "CET", multiple = FALSE),
    hr(),
    textOutput("status"),
    div(style = "color: #ff6b6b; font-size: 0.8em;", textOutput("Messages_Resu"))
  ),
  
  # --- En-tête personnalisé avec Logo (à mettre dans le corps de la page) ---
  div(
    style = "display: flex; justify-content: space-between; align-items: center; padding: 10px 20px; background-color: #0f0f0f; border-bottom: 1px solid #333; margin-bottom: 20px;",
    div(
      style = "font-size: 24px; font-weight: bold; color: #fff;",
      "BirdNET-ResChecker - PHENOLOGY"
    ),
    div(
      tags$img(
        src = "logo.png", 
        height = "70px", 
        style = "margin-right: 10px; border-radius: 12px; padding: 6px 10px; background-color: white;"
      ) 
    )
  ),
  
  navset_card_underline(
    nav_panel("Phenology",
              layout_columns(
                col_widths = c(4, 8),
                # Card 7: Pheno Params
                card(full_screen = TRUE, card_header("Phenology Parameters"), 
                     selectizeInput(inputId = "species_pheno", label = "Species List",
                                    choices = NULL, selected = "All species", multiple = TRUE,
                                    options = list(placeholder = 'Select species or "All Species"', create = FALSE)),
                     numericInput("Confid_Pheno", "Confidence Threshold", min = 0.01, max = 1, value = 0.01),
                     layout_columns(col_widths = c(6, 6),
                                    numericInput("Lat", "Latitude", min = -90, max = 90, value = 46.97),
                                    numericInput("Lon", "Longitude", min = -180, max = 180, value = 6.97)),
                     sliderInput("Unit", "Aggregation Interval (min)", min = 1, max = 60, value = 15),
                     checkboxInput("Noctu_plot", label = "Nocturnal Plot?", value = FALSE),
                     actionButton("start_pheno", "Generate Plot", class = "btn-primary", width = "100%")),
                # Card 8: Pheno Plot
                card(full_screen = TRUE, card_header("Phenology Graph"), plotOutput("pheno_plot"))
              )),
    
    nav_panel("Reference",
              layout_columns(
                col_widths = c(2, 8, 2),
                textOutput("description"),
                textOutput("contributions"),
                textOutput("license")
              ))
  ),
  
  # CSS Personnalisé pour adapter les composants spécifiques (inputs, tableaux)
  tags$head(
    tags$style(HTML("
      /* Fond général et texte */
      body { background-color: #2b2b2b; color: #e0e0e0; }
      
      /* Inputs (Selectize, Text, Date) */
      .form-control, .selectize-input, .selectize-control.multi .selectize-input > div {
        background-color: #444444 !important;
        color: #ffffff !important;
        border: 1px solid #555555 !important;
      }
      .selectize-dropdown, .selectize-dropdown-content {
        background-color: #444444 !important;
        color: #ffffff !important;
      }
      .selectize-dropdown .option { color: #ffffff !important; }
      .selectize-dropdown .active { background-color: #555555 !important; color: #ffffff !important; }
      
      /* Labels */
      .control-label, .bslib-sidebar-input label { color: #ffffff !important; font-weight: 600; }
      
      /* Tableaux DT */
      .dataTables_wrapper .dataTables_length, 
      .dataTables_wrapper .dataTables_filter, 
      .dataTables_wrapper .dataTables_info, 
      .dataTables_wrapper .dataTables_processing, 
      .dataTables_wrapper .dataTables_paginate { color: #e0e0e0 !important; }
      
      .dataTable tbody tr { color: #e0e0e0; }
      .dataTable tbody tr:hover { background-color: #3a3a3a; }
      
      /* Cartes et Panneaux */
      .card { background-color: #333333; border: 1px solid #444444; color: #ffffff; }
      .card-header { background-color: #3a3a3a; border-bottom: 1px solid #444444; color: #ffffff; }
      
      /* Scrollbars */
      ::-webkit-scrollbar { width: 10px; }
      ::-webkit-scrollbar-track { background: #2b2b2b; }
      ::-webkit-scrollbar-thumb { background: #555555; border-radius: 5px; }
      ::-webkit-scrollbar-thumb:hover { background: #777777; }
    "))
  )
)

# --- 5. LOGIQUE SERVEUR (SERVER) ----------------------------------------------
server <- function(input, output, session) {
  
  # --- Textes statiques ---
  output$description <- renderText("This application analyzes BirdNET detection results. Select your folders on the left, then explore detections, spectrograms, and phenology graphs.")
  output$contributions <- renderText("Contributions: Amandine Serrurier, Jean-Nicolas Pradervand & Christophe Sahli\nSwiss Ornithological Institute")
  output$license <- renderText("MIT License © 2025 Christophe Sahli")
  
  # --- 5.1 Gestion des Dossiers (ShinyFiles) ---
  # fileInput(
  #   "dir1_files",
  #   "Choose files – folder 1",
  #   multiple = TRUE
  # )
  # 
  # fileInput(
  #   "dir2_files",
  #   "Choose files – folder 2",
  #   multiple = TRUE
  # )
  # shinyDirChoose(input, "dir1", roots = c(home = roots_home))
  # shinyDirChoose(input, "dir2", roots = c(home = roots_home))
   
  # observeEvent(input$dir1, {
  #   tryCatch({
  #     dir1_path <- parseDirPath(roots = c(home = roots_home), input$dir1)
  #     if (is.character(dir1_path) && length(dir1_path) > 0 && dir1_path != "") {
  #       dir1(dir1_path)
  #     }
  #   }, error = function(e) {
  #     showNotification("Invalid path for results folder.", type = "warning")
  #   })
  # })
  # 
   # 
  
  # dir2 <- reactiveVal(NULL)
  
  # observeEvent(input$dir1, {
  #   req(input$dir1)
  #   
  # })
   dir1 <- reactiveVal(NULL)
   observeEvent(input$dir1, {
     
     files <- input$dir1
     req(files)
     
     if (identical(attr(files, "upload_status"), "completed")) {
       
       dir1(attr(files, "upload_dir"))
       
     }
   })
  
  # output$dir1_path <- renderText({
  #   if(!is.null(dir1()) && dir1() != "") paste("Results:", basename(dir1())) else "No folder selected"
  # })
  # output$dir2_path <- renderText({
  #   if(!is.null(dir2()) && dir2() != "") paste("Audio:", basename(dir2())) else "No folder selected"
  # })
  # 
  # --- 5.2 Chargement des Données ---
  DT_reac <- reactiveVal(NULL)
  Sp_List <- reactiveVal(NULL)
  selected_list_species <- reactiveVal(NULL)
  complete_audio_path <- reactiveVal(NULL)
  Date_DT_min <- reactiveVal(NULL)
  Date_DT_max <- reactiveVal(NULL)
  
  # observeEvent(input$start, {
  #   req(input$dir1)
  #   DT_reac(NULL)
  #   showNotification("Loading selection tables...", closeButton = FALSE, duration = 2)
  #   
  #   if (!exists("load_selection_tables")) {
  #     showNotification("Error: Function 'load_selection_tables' not found.", type = "error")
  #     return()
  #   }
  observeEvent(input$start, {
    
    DT_reac(NULL)
    
    showNotification(
      "Loading selection tables...",
      closeButton = FALSE,
      duration = 2
    )
    
    tryCatch({
      
      if (isTRUE(input$Compiled_F)) {
        
        # ---- COMPILED FORMAT ----
        req(input$compiled_file)
        
        path1 <- input$compiled_file$datapath
        
        # Lecture directe du fichier
        # DT <- data.table::fread(path1)
        # names(DT) <- make.names(names(DT))
       DT <- load_selection_tables_ONLINE(
          dir1 = path1,
          dir2 = NULL,
          compiled = TRUE,
          utc_tz = input$UTC_choice
        )
        print(names(DT))
        print(str(DT))
        
      } else {
        
        # ---- NON COMPILED FORMAT ----
        req(dir1())
        
        path1 <- dir1()
        
        # On garde exactement ton fonctionnement précédent
        DT <- load_selection_tables_ONLINE(
          dir1 = path1,
          dir2 = NULL,
          compiled = FALSE,
          utc_tz = input$UTC_choice
        )
        print(names(DT))
        print(str(DT))
        
      }
      
      # ---- À partir d'ici, DT existe dans les deux cas ----
      
      if (is.null(DT) || nrow(DT) == 0) {
        showNotification(
          "No data found or empty tables.",
          type = "error"
        )
        output$status <- renderText("No data loaded.")
        return()
      }
      
      # complete_audio_path(DT$Begin_Path)
      print("names after treatment")
      print(names(DT))
      
      Date_DT_min(min(DT$Date, na.rm = TRUE))
      Date_DT_max(max(DT$Date, na.rm = TRUE))
      
      updateDateRangeInput(
        session,
        "filter_dates",
        start = Date_DT_min(),
        end = Date_DT_max(),
        min = Date_DT_min(),
        max = Date_DT_max()
      )
      
      Sp_List_DT <- sort(unique(DT$Common.Name))
      Sp_List(Sp_List_DT)
      
      species_list_pheno <- c(Sp_List_DT, "All species")
      
      updateSelectizeInput(
        session,
        "species_pheno",
        choices = species_list_pheno,
        selected = "All species",
        server = TRUE
      )
      
      DT_reac(DT)
      
      output$status <- renderText(
        paste(
          "Loaded:",
          nrow(DT),
          "detections |",
          length(Sp_List_DT),
          "species"
        )
      )
      
    }, error = function(e) {
      
      showNotification(
        paste("Error loading data:", e$message),
        type = "error"
      )
      
      output$status <- renderText("Error loading data.")
      
    })
    
  })
#   observeEvent(input$start, {
#     
#     if (isTRUE(input$Compiled_F)) {
#       
#       req(input$compiled_file)
#       
#       path1 <- input$compiled_file$datapath
#       
#     } else {
#       
#       req(dir1())
#       
#       path1 <- dir1()
#       
#     }
#     
#     DT_reac(NULL)
#     
#     showNotification(
#       "Loading selection tables...",
#       closeButton = FALSE,
#       duration = 2
#     )
#     
#     if (!exists("load_selection_tables")) {
#       showNotification(
#         "Error: Function 'load_selection_tables' not found.",
#         type = "error"
#       )
#       return()
#     }
#     
#     tryCatch({
#       
#       DT <- load_selection_tables(
#         dir1 = path1,
#         dir2 = dir2(),
#         compiled = input$Compiled_F,
#         utc_tz = input$UTC_choice
#       )
#       
#     tryCatch({
#       DT <- load_selection_tables(
#         dir1 = dir1(), dir2 = dir2(), compiled = input$Compiled_F, utc_tz = input$UTC_choice
#       )
# 
#    
# 
#       if (is.null(DT) || nrow(DT) == 0) {
#         showNotification("No data found or empty tables.", type = "error")
#         output$status <- renderText("No data loaded.")
#         return()
#       }
#       
#       complete_audio_path(DT$Begin_Path)
#       Date_DT_min(min(DT$Date, na.rm = TRUE))
#       Date_DT_max(max(DT$Date, na.rm = TRUE))
#       
#       updateDateRangeInput(session, "filter_dates", 
#                            start = Date_DT_min(), end = Date_DT_max(),
#                            min = Date_DT_min(), max = Date_DT_max())
#       
#       Sp_List_DT <- sort(unique(DT$Common.Name))
#       Sp_List(Sp_List_DT)
#       
#       species_list_pheno <- c(Sp_List_DT, "All species")
#       updateSelectizeInput(session, "species_pheno", choices = species_list_pheno, 
#                            selected = "All species", server = TRUE)
#       
#       DT_reac(DT)
#       output$status <- renderText(paste("Loaded:", nrow(DT), "detections |", length(Sp_List_DT), "species"))
#       
#     }, error = function(e) {
#       showNotification(paste("Error loading data:", e$message), type = "error")
#       output$status <- renderText("Error loading data.")
#     })
#   })
#   
  # --- 5.3 Filtrage Réactif ---
  DT_filtered <- reactiveVal(NULL)
  last_species_list <- reactiveVal(NULL)

  filtered_confid <- reactive(input$Confid) %>% debounce(300)
  filtered_dates <- reactive(input$filter_dates) %>% debounce(300)

  observe({
    req(DT_reac())
    req(filtered_confid(), filtered_dates())

    DT_filter <- copy(DT_reac())

    conf_val <- filtered_confid()
    DT_filter <- DT_filter[Confidence >= conf_val[1] & Confidence <= conf_val[2]]

    dates_val <- filtered_dates()
    if (!is.null(dates_val) && length(dates_val) == 2 && !any(is.na(dates_val))) {
      DT_filter <- DT_filter[Date >= dates_val[1] & Date <= dates_val[2]]
    }

    current_species <- sort(unique(DT_filter$Common.Name))
    selected_list_species(current_species)

    current_sel <- input$species
    species_valid <- !is.null(current_sel) && current_sel != "" && current_sel %in% current_species

    if (species_valid) {
      DT_filter <- DT_filter[Common.Name %in% current_sel]
    }

    old_list <- last_species_list()
    needs_update <- FALSE

    if (is.null(old_list) || length(old_list) != length(current_species)) {
      needs_update <- TRUE
    } else if (!identical(old_list, current_species)) {
      needs_update <- TRUE
    } else if (!species_valid && !is.null(current_sel)) {
      needs_update <- TRUE
    }

    if (needs_update) {
      updateSelectizeInput(
        session,
        "species",
        choices = current_species,
        selected = if(species_valid) current_sel else NULL,
        server = TRUE
      )
      last_species_list(current_species)
    }

    DT_filtered(DT_filter)
  })

  # --- 5.8 Phénologie ---
  pheno_data <- eventReactive(input$start_pheno, {
    req(DT_reac())
    if (!exists("pheno_matrix")) {
      showNotification("Function 'pheno_matrix' missing.", type = "error")
      return(NULL)
    }

    dt_pheno <- DT_reac()
    xlim_plot <- c(as.Date(min(dt_pheno$Date, na.rm = TRUE)),
                   as.Date(max(dt_pheno$Date, na.rm = TRUE)))

    sp_to_plot <- input$species_pheno
    if (is.null(sp_to_plot) || length(sp_to_plot) == 0 || "All species" %in% sp_to_plot) {
      sp_to_plot <- unique(dt_pheno$Common.Name)
    }

    tryCatch({
      pheno_matrix(
        Voc = dt_pheno, SP = sp_to_plot, Unit = input$Unit,
        Confidence1 = input$Confid_Pheno, sunrise = TRUE,
        LAT = input$Lat, LONG = input$Lon, UTC = input$UTC_choice,
        xlim_plot = xlim_plot, nocturnal = input$Noctu_plot
      )
    }, error = function(e) {
      showNotification(paste("Pheno Error:", e$message), type = "error")
      return(NULL)
    })
  })

  output$pheno_plot <- renderPlot({
    data <- pheno_data()
    if (is.null(data)) {
      plot(1, type="n", axes=FALSE, xlab="", ylab="",
           main = "No data to plot\n(Lower confidence threshold or select species)")
      text(1, 1, "Check parameters", cex=1.2, col="red")
    } else {
      if (inherits(data, "gg")) print(data) else plot(data)
    }
  })
}
# --- 6. LANCEMENT ---
shinyApp(ui = ui, server = server)