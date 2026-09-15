# ==============================================================================
# BirdNET-ResChecker - PHENOLOGY
# (High Memory & Config)
# ==============================================================================

# --- 1. CHARGEMENT DES LIBRAIRIES ---------------------------------------------
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
  library(pacman)
}

pacman::p_load(
  "shiny", "grid", "bslib", "stringr", "data.table", "shinyFiles", "tidyr", "janitor",
  "stringi", "fs", "parallel", "ggplot2", "mgcv", "gratia", "hms", "suncalc", "scales",
  "lubridate", "DT", "future.apply", "here"
)


# --- 2. CHARGEMENT DES FONCTIONS PERSONNALISÉES -------------------------------
# NOTE CRITIQUE POUR GROS FICHIERS :
# La fonction 'load_selection_tables.R' doit utiliser data.table::fread() 
# et NON read.csv() ou read.table() pour être efficace sur des fichiers > 500Mo.
if (dir.exists("Common_functions")) {

  source("Common_functions/get_roots.R")
  source("Common_functions/pheno_matrix.R")
  source("Common_functions/clean_text.R")
  source("Common_functions/export_selected_audio.R")
  source("Common_functions/prepare_for_spectro.R")
  source("Common_functions/load_selection_tables.R")
  source("Common_functions/BirdNET_data_parser.R")
} else {
  warning("Folder 'Common_functions' not found. Simple mode activated.")
  get_roots <- function() return("~") 
  load_selection_tables <- function(...) return(NULL)
  prepare_for_spectro <- function(...) return(list(NULL, 22050, NULL, 16))
  pheno_matrix <- function(...) return(NULL)
  export_selected_audio <- function(wav, ...) return(wav)
  clean_text <- function(x) return(gsub("[^a-zA-Z0-9]", "_", x))
}

# --- 3. CONFIGURATION GLOBALE -------------------------------------------------
roots_home <- get_roots()
utc_timezones <- OlsonNames()

# Augmenter la limite de mémoire pour R (si le système le permet)
# Utile pour les sessions locales. Sur serveur, dépend de la config RAM allouée.
if (!grepl("shiny-server", Sys.info()["user"])) {
  memory.limit(size = 50000) # Tente d'augmenter la limite à 50GB sur Windows
}

dark_theme <- bs_theme(
  version = 5,
  bootswatch = "darkly",
  bg = "#1a1a1a",
  fg = "#e0e0e0",
  primary = "#375a7f",
  secondary = "#444444"
)

ui <- page_sidebar(
  theme = dark_theme,
  sidebar = sidebar(
    title = "Data Source",
    shinyDirButton("dir1", "Select BirdNET Results Folder", title = "Choose Results Folder"),
    checkboxInput("Compiled_F", label = "Results in compiled format?", value = FALSE),
    textOutput("dir1_path"),
    actionButton("start", "Start App", class = "btn-primary", width = "100%"),
    hr(),
    selectInput(inputId = "device_tz", label = "Timezone of set device", 
                choices = utc_timezones, selected = Sys.timezone(), multiple = FALSE),
    selectInput(inputId = "deployment_tz", label = "Timezone of deployement", 
                choices = utc_timezones, selected = Sys.timezone(), multiple = FALSE),
    actionButton("restart_timezone", "Update Timezones", class = "btn-primary", width = "100%"),
    hr(),
    uiOutput("recorder_ui"),
    hr(),
    textOutput("status"),
    div(style = "color: #ff6b6b; font-size: 0.8em;", textOutput("Messages_Resu"))
  ),
  
  div(
    style = "display: flex; justify-content: space-between; align-items: center; padding: 10px 20px; background-color: #0f0f0f; border-bottom: 1px solid #333; margin-bottom: 20px;",
    div(style = "font-size: 24px; font-weight: bold; color: #fff;", "Bioacoustic tools - PHENOLOGY"),
    div(tags$img(src = "logo.png", height = "70px", style = "margin-right: 10px; border-radius: 12px; padding: 6px 10px; background-color: white;"))
  ),
  
  navset_card_underline(
    nav_panel("Phenology",
              layout_columns(
                col_widths = c(4, 8),
                card(full_screen = TRUE, card_header("Phenology Parameters"),
                     uiOutput("species_pheno_ui"),
                     numericInput("Confid_Pheno", "Confidence Threshold", min = 0.01, max = 1, value = 0.01),
                     layout_columns(col_widths = c(6, 6),
                                    numericInput("Lat", "Latitude", min = -90, max = 90, value = 46.97),
                                    numericInput("Lon", "Longitude", min = -180, max = 180, value = 6.97)),
                     sliderInput("Unit", "Aggregation Interval (min)", min = 1, max = 60, value = 15),
                     checkboxInput("Noctu_plot", label = "Nocturnal Plot?", value = FALSE),
                     actionButton("start_pheno", "Generate Plot", class = "btn-primary", width = "100%")),
                card(full_screen = TRUE, card_header("Phenology Graph"), plotOutput("pheno_plot"))
              )),
    nav_panel("Reference",
              layout_columns(col_widths = c(2, 8, 2),
                             textOutput("description"),
                             textOutput("contributions"),
                             textOutput("license")))
  ),
  
  tags$head(
    tags$style(HTML("
      body { background-color: #2b2b2b; color: #e0e0e0; }
      .form-control, .selectize-input, .selectize-control.multi .selectize-input > div {
        background-color: #444444 !important; color: #ffffff !important; border: 1px solid #555555 !important;
      }
      .selectize-dropdown, .selectize-dropdown-content {
        background-color: #444444 !important; color: #ffffff !important;
      }
      .selectize-dropdown .option { color: #ffffff !important; }
      .selectize-dropdown .active { background-color: #555555 !important; color: #ffffff !important; }
      .control-label, .bslib-sidebar-input label { color: #ffffff !important; font-weight: 600; }
      .dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter, 
      .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_processing, 
      .dataTables_wrapper .dataTables_paginate { color: #e0e0e0 !important; }
      .dataTable tbody tr { color: #e0e0e0; }
      .dataTable tbody tr:hover { background-color: #3a3a3a; }
      .card { background-color: #333333; border: 1px solid #444444; color: #ffffff; }
      .card-header { background-color: #3a3a3a; border-bottom: 1px solid #444444; color: #ffffff; }
      ::-webkit-scrollbar { width: 10px; }
      ::-webkit-scrollbar-track { background: #2b2b2b; }
      ::-webkit-scrollbar-thumb { background: #555555; border-radius: 5px; }
      ::-webkit-scrollbar-thumb:hover { background: #777777; }
    "))
  )
)

# --- 5. LOGIQUE SERVEUR (SERVER) ----------------------------------------------
server <- function(input, output, session) {
  
  output$description <- renderText("This application analyzes BirdNET detection results and produces phenology graphs.")
  output$contributions <- renderText("Contributions: Amandine Serrurier, Jean-Nicolas Pradervand & Christophe Sahli\nSwiss Ornithological Institute")
  output$license <- renderText("MIT License © 2026 Swiss Ornithological Institute")
  
  shinyDirChoose(input, "dir1", roots = c(home = roots_home))
  shinyDirChoose(input, "dir2", roots = c(home = roots_home))
  
  dir1 <- reactiveVal(NULL)
  dir2 <- reactiveVal(NULL)
  
  observeEvent(input$dir1, {
    tryCatch({
      dir1_path <- parseDirPath(roots = c(home = roots_home), input$dir1)
      if (is.character(dir1_path) && length(dir1_path) > 0 && dir1_path != "") {
        dir1(dir1_path)
      }
    }, error = function(e) {
      showNotification("Invalid path for results folder.", type = "warning")
    })
  })
  
  observeEvent(input$dir2, {
    tryCatch({
      dir2_path <- parseDirPath(roots = c(home = roots_home), input$dir2)
      if (is.character(dir2_path) && length(dir2_path) > 0 && dir2_path != "") {
        dir2(dir2_path)
      }
    }, error = function(e) {
      showNotification("Invalid path for audio folder.", type = "warning")
    })
  })
  
  output$dir1_path <- renderText({
    if(!is.null(dir1()) && dir1() != "") paste("Results:", basename(dir1())) else "No folder selected"
  })
  output$dir2_path <- renderText({
    if(!is.null(dir2()) && dir2() != "") paste("Audio:", basename(dir2())) else "No folder selected"
  })
  
  DT_reac <- reactiveVal(NULL)
  Sp_List <- reactiveVal(NULL)
  selected_list_species <- reactiveVal(NULL)
  complete_audio_path <- reactiveVal(NULL)
  Date_DT_min <- reactiveVal(NULL)
  Date_DT_max <- reactiveVal(NULL)
  
  observeEvent(input$start, {
    req(input$dir1)
    DT_reac(NULL)
    
    # Utilisation de withProgress pour les gros chargements
    withProgress(message = "Loading large dataset...", detail = "Parsing selection tables, please wait.", value = 0, {
      
      if (!exists("load_selection_tables")) {
        showNotification("Error: Function 'load_selection_tables' not found.", type = "error")
        return()
      }
      
      tryCatch({
        # Simulation de progression si le chargement est long
        incProgress(0.2, detail = "Initializing...")
        
        DT <- load_selection_tables(
          dir1 = dir1(), 
          dir2 = dir2(), 
          compiled = input$Compiled_F, 
          device_tz = input$device_tz,      # Ex: "America/Anchorage" (choisi par l'utilisateur)
          deployment_tz = input$deployment_tz # Ex: "Europe/Zurich" (choisi par l'utilisateur)
        )
        
        incProgress(0.5, detail = "Finalizing data structure...")
        
        if (is.null(DT) || nrow(DT) == 0) {
          showNotification("No data found or empty tables.", type = "error")
          output$status <- renderText("No data loaded.")
          return()
        }
        
        # Nettoyage mémoire immédiat après gros chargement
        gc(verbose = FALSE)
        incProgress(0.3, detail = "Memory optimization...")
        
        complete_audio_path(DT$Begin_Path)
        Date_DT_min(min(DT$Date, na.rm = TRUE))
        Date_DT_max(max(DT$Date, na.rm = TRUE))
        
        updateDateRangeInput(session, "filter_dates", 
                             start = Date_DT_min(), end = Date_DT_max(),
                             min = Date_DT_min(), max = Date_DT_max())
        
        Sp_List_DT <- sort(unique(DT$Common.Name))
        
        setkey(DT, Common.Name, Date, Confidence) # créer des keys pour un tri très rapide
        DT_reac(DT)
        output$status <- renderText(paste("Loaded:", nrow(DT), "detections |", length(Sp_List_DT), "species"))
        
      }, error = function(e) {
        msg <- e$message
        if (grepl("cannot allocate", msg, ignore.case = TRUE)) {
          showNotification("ERROR: Not enough RAM to load this file. Please increase server memory or split data.", type = "error", duration = 10)
          output$status <- renderText("Error: Memory Limit Reached")
        } else {
          showNotification(paste("Error loading data:", msg), type = "error")
          output$status <- renderText("Error loading data.")
        }
      })
    })
  })
  
  
  
  observeEvent(input$restart_timezone, {
    req(DT_reac())
    req(input$device_tz)
    req(input$deployment_tz)
    
    # 1. CRÉER UNE COPIE EXPLICITE POUR ÉVITER LA MODIFICATION DIRECTE
    DT_new <- data.table::copy(DT_reac())
    
    showNotification("Recalculating timezones...\nAnother message will appear when it's done.", type = "message", duration = 4)
    
    tryCatch({
      
      # Création de DateTime_Display dans le fuseau de l'appareil (ex: Alaska)
      DT_new[, DateTime_Display := as.POSIXct(date_strings, format = "%Y%m%d_%H%M%S", tz = input$device_tz)]
      #attr(DT_new$DateTime_Display, "tzone") <- device_tz # Verrouillage
      
      # Colonnes dérivées "Display" (pour référence/logs)
      DT_new[, `:=`(
        Date_Display = as.Date(DateTime_Display),
        Hour_Display = hour(DateTime_Display),
        Min_Display = minute(DateTime_Display),
        Time_Display = sprintf("%02d:%02d", hour(DateTime_Display), minute(DateTime_Display))
      )]
      
      # Gestion des Dates - CONVERSION VERS L'HEURE "RÉELLE" (DEPLOYMENT)
      if (!is.null(input$deployment_tz) && input$deployment_tz != input$device_tz) {
        # Conversion via UTC pour garantir la justesse de l'instant
        # 1. On passe en UTC (instant universel)
        utc_strings <- format(DT_new$DateTime_Display, tz = "UTC", usetz = FALSE)
        
        # 2. On recrée l'objet dans le fuseau de déploiement (ex: France)
        DT_new[, DateTime_Real := as.POSIXct(utc_strings, format = "%Y-%m-%d %H:%M:%S", tz = input$deployment_tz)]
        attr(DT_new$DateTime_Real, "tzone") <- input$deployment_tz
        
        message("Conversion des fuseaux horaires effectuée.")
      } else {
        # Si mêmes fuseaux, on duplique
        DT_new[, DateTime_Real := DateTime_Display]
        message("Pas de conversion de fuseau (identiques).")
      }
      
      
      # Calcul du segment temporel absolu (Start_segment)
      # On utilise l'heure RÉELLE pour les calculs scientifiques (soleil, etc.)
      if (input$Compiled_F) {
        DT_new[, time_offset := File.Offset..s.]
      } else {
        if ("Begin.Time..s." %in% names(DT_new)) {
          DT_new[, time_offset := Begin.Time..s.]
        } else {
          DT_new[, time_offset := 0]
        }
      }
      
      # Calcul de Start_segment basé sur DateTime_Real
      DT_new[, Start_segment := DateTime_Real + time_offset]
      
      # Calcul Stop_segment
      if (input$Compiled_F) {
        if ("End.Time..s." %in% names(DT_new)) {
          DT_new[, Stop_segment := Start_segment + (End.Time..s. - File.Offset..s.)]
        } else {
          DT_new[, Stop_segment := Start_segment]
        }
      } else {
        if ("End.Time..s." %in% names(DT_new)) {
          DT_new[, Stop_segment := DateTime_Real + End.Time..s.]
        } else {
          DT_new[, Stop_segment := Start_segment]
        }
      }
      
      DT_new[, `:=`(
        Date = as.Date(Start_segment),
        Hour = hour(Start_segment),
        Min = minute(Start_segment),
        Time = sprintf("%02d:%02d", hour(Start_segment), minute(Start_segment)),
        Hour_Decimal_Real = hour(Start_segment) + (minute(Start_segment) / 60)
      )]
      
      DT_reac(DT_new)
      
      showNotification(paste("Timezone updated:", input$device_tz, "->", input$deployment_tz), type = "message", duration = 2)
      
    }, error = function(e) {
      showNotification(paste("Error updating timezone:", e$message), type = "error")
      # En cas d'erreur, DT_reac n'est pas modifié, on garde les données précédentes
    })
  })
  
  
  

  # ---  UI DYNAMIQUE ---
  output$recorder_ui <- renderUI({
    req(DT_reac())
    df <- DT_reac()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    recs <- sort(unique(na.omit(df$recorder)))
    if(length(recs) == 0) return(NULL)
    selectInput("selected_recorders", "Select recorder(s):", choices = recs, selected = recs, multiple = TRUE)
  })
  
  output$species_pheno_ui <- renderUI({
    req(DT_reac())
    df <- DT_reac()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    species <- sort(unique(na.omit(df$Common.Name)))
    if(length(species) == 0) return(NULL)
    selectInput("species_pheno", "Select species", choices = c("All species", species), selected = "All species", multiple = TRUE)
  })
  
  
  # --- 5.8 Phénologie ---
  pheno_data <- eventReactive(input$start_pheno, {
    req(DT_reac())
    print(head(DT_reac())) # debeug
    if (!exists("pheno_matrix")) {
      showNotification("Function 'pheno_matrix' missing.", type = "error")
      return(NULL)
    }
    
    dt_pheno <- DT_reac()
    
    print("dt_pheno") # debeug
    print(head(dt_pheno)) # debeug
    
    xlim_plot <- c(as.Date(min(dt_pheno$Date, na.rm = TRUE)), 
                   as.Date(max(dt_pheno$Date, na.rm = TRUE)))
    
    print("xlim_plot") # debeug
    print(xlim_plot) # debeug
    
    conf_val <- input$Confid_Pheno

    print("conf_val") # debeug
    print(conf_val) # debeug
    
    # 2. Filtrage Intermédiaire
    # C'est sur ce jeu de données qu'on calculera la liste des espèces et recorders disponibles
    DT_intermediate <- dt_pheno[
      Confidence >= conf_val
    ]
    
    # if (!is.null(dates_val) && length(dates_val) == 2 && !any(is.na(dates_val))) {
    #   DT_intermediate <- DT_intermediate[Date >= dates_val[1] & Date <= dates_val[2]]
    # }
    
    # 3. Filtrage par Enregistreurs
    # On vérifie si l'input existe et n'est pas vide
    if (!is.null(input$selected_recorders) && length(input$selected_recorders) > 0) {
      DT_intermediate <- DT_intermediate[recorder %in% input$selected_recorders]
    }
    
    
    print("DT_intermediate") # debeug
    print(head(DT_intermediate)) # debeug
    
    
    sp_to_plot <- input$species_pheno

    print("sp_to_plot") # debeug
    print(sp_to_plot) # debeug
    
    
      # FILTRAGE PRIORITAIRE PAR ESPÈCE
      if (!"All species" %in% sp_to_plot) {
        Voc <- DT_intermediate[Common.Name %in% sp_to_plot]
      } else {
        Voc <- DT_intermediate
      }
    
    
    print("Voc avant matrix") # debeug
    print(head(Voc)) # debeug
    
    
    tryCatch({
      pheno_matrix(
        Voc = Voc, SP = sp_to_plot, Unit = input$Unit, 
        Confidence1 = input$Confid_Pheno, sunrise = TRUE, 
        LAT = input$Lat, LONG = input$Lon, TimeZone = input$deployment_tz,
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

# --- LANCEMENT ---
shinyApp(ui = ui, server = server)