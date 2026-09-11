

# ==============================================================================
# BirdNET-ResChecker - SPECTRO
# Version Optimisée & Corrigée
# ==============================================================================

# --- 1. CHARGEMENT DES LIBRAIRIES ---------------------------------------------
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
  library(pacman)
}

# Installation et chargement des paquets requis
pacman::p_load(
  "shiny", "grid", "bslib", "stringr", "data.table", "shinyFiles", "tidyr", 
  "stringi", "fs", "parallel", "ggplot2", "mgcv", "gratia", "hms", "suncalc",
  "lubridate", "DT", "seewave", "tuneR", "base64enc", "av", "audio", "future.apply", "signal"
)

# --- 2. CHARGEMENT DES FONCTIONS PERSONNALISÉES -------------------------------
# Vérifie si le dossier existe avant de sourcer les fichiers
if (dir.exists("Common_functions")) {
  source("Common_functions/get_roots.R")
  source("Common_functions/pheno_matrix.R")
  source("Common_functions/clean_text.R")
  source("Common_functions/export_selected_audio.R")
  source("Common_functions/prepare_for_spectro.R")
  source("Common_functions/load_selection_tables.R")
} else {
  # Fallback pour éviter le crash si le dossier est manquant (mode debug)
  warning("Dossier 'Common_functions' non trouvé. Assurez-vous qu'il est dans le répertoire de travail.")
  # Définition de fonctions factices pour que l'UI se charge (à retirer en prod)
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

# Messages d'erreur standardisés
Files_message_1 <- "-- There are no results files (~.BirdNET.selection.table.txt) in the results folder. --"
Files_message_2 <- "-- There are no audio files in the audio folder. --"
Files_message_3 <- "-- Audio files are missing from the audio folder. --"
# --- 4. INTERFACE UTILISATEUR (UI) --------------------------------------------

# Définition du thème sombre avec bslib
dark_theme <- bs_theme(
  version = 5,
  bootswatch = "darkly", # Thème sombre prédéfini
  bg = "#1a1a1a",        # Fond principal
  fg = "#e0e0e0",        # Couleur du texte
  primary = "#375a7f",   # Couleur primaire (bleu sombre)
  secondary = "#444444"  # Couleur secondaire
)

ui <- page_fluid( # On passe à page_fluid pour mieux contrôler l'en-tête
  theme = dark_theme,
  
  # --- En-tête personnalisé avec Logo ---
  div(
    style = "display: flex; justify-content: space-between; align-items: center; padding: 10px 20px; background-color: #0f0f0f; border-bottom: 1px solid #333; margin-bottom: 20px;",
    div(
      style = "font-size: 24px; font-weight: bold; color: #fff;",
      "BirdNET-ResChecker - SPECTRO"
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
  
  # --- Sidebar et Contenu ---
  layout_sidebar(
    sidebar = sidebar(
      title = "Data Source",
      shinyDirButton("dir1", "Select BirdNET Results Folder", title = "Choose Results Folder"),
      checkboxInput("Compiled_F", label = "Results in compiled format?", value = FALSE),
      textOutput("dir1_path"),
      hr(),
      shinyDirButton("dir2", "Select Audio Folder", title = "Choose Audio Folder"),
      textOutput("dir2_path"),
      hr(),
      selectInput(inputId = "device_tz", label = "Timezone of set device", 
                  choices = utc_timezones, selected = Sys.timezone(), multiple = FALSE),
      selectInput(inputId = "deployment_tz", label = "Timezone of deployement", 
                  choices = utc_timezones, selected = Sys.timezone(), multiple = FALSE),
      hr(),
      actionButton("start", "Start App", class = "btn-primary", width = "100%"),
      hr(),
      textOutput("status"),
      div(style = "color: #ff6b6b; font-size: 0.8em;", textOutput("Messages_Resu")),
      div(style = "color: #ffa500; font-size: 0.8em;", textOutput("Messages_Audio"))
    ),
    
    # Contenu principal
    navset_card_underline(
      id = "tabs", # ID requis pour la stabilité
      
      nav_panel("Spectrogram & Export", 
                layout_columns(
                  col_widths = c(8, 4),
                  card(full_screen = TRUE, card_header("Spectrogram & Audio"), 
                       uiOutput("Spectro_error"), # Ajout pour les messages d'erreur
                       plotOutput("Spectro"),
                       uiOutput("audio_player")),
                  card(full_screen = TRUE, card_header("Parameters"), 
                       sliderInput("Highpass", "High-Pass Filter (Hz)", min = 0, max = 8000, value = 200),
                       sliderInput("Freq_Spectro", "Frequency Range (Y-axis, Hz)", min = 0, max = 12000, value = c(0, 10000), dragRange = TRUE))
                  ),
                layout_columns(
                  col_widths = c(3, 9),
                  card(full_screen = TRUE, card_header("Filters"), 
                       sliderInput("Confid", "Confidence Interval", min = 0.01, max = 1, value = c(0.01, 1), dragRange = TRUE),
                       selectizeInput(inputId = "species", label = "Filter by Species",
                                      choices = NULL, multiple = FALSE,
                                      options = list(placeholder = 'Select a species', create = FALSE)),
                       dateRangeInput(inputId = "filter_dates", label = "Date Range",
                                      start = as.Date("1900-01-01"), end = as.Date(Sys.Date())),
                       sliderInput("Hour_Range", "Hour range", min = 0, max = 23, value = c(0, 23), dragRange = TRUE),
                       uiOutput("recorder_ui"),
                       textOutput("selected_audio")),
                  card(full_screen = TRUE, card_header("Actions & Data"), 
                       div(style = "margin-bottom: 15px;",
                           downloadButton('download', "Export CSV", class = "btn-success"),
                           actionButton("open_audio", "Open Audio File", class = "btn-info"),
                           downloadButton('download_audio', "Export WAV", class = "btn-secondary"),
                           downloadButton('download_audio_mp3', "Export MP3", class = "btn-secondary")
                       ),
                       DTOutput("summary"))
                )
                ),
      
      nav_panel("Reference",
                layout_columns(
                  col_widths = c(2, 8, 2),
                  textOutput("description"),
                  textOutput("contributions"),
                  textOutput("license")
                ))
    )
  ),
  
  # --- CSS Personnalisé pour forcer le mode sombre partout ---
  tags$head(
    tags$style(HTML("
      /* Fond général et texte */
      body { background-color: #1a1a1a; color: #e0e0e0; }
      
      /* Cartes et panneaux */
      .card, .bslib-card { background-color: #2c2c2c !important; border-color: #444 !important; }
      .card-header { background-color: #333 !important; color: #fff !important; border-bottom: 1px solid #555; }
      
      /* Tableaux DT (DataTables) - Important pour le contraste */
      .dataTables_wrapper .dataTables_length, 
      .dataTables_wrapper .dataTables_filter, 
      .dataTables_wrapper .dataTables_info, 
      .dataTables_wrapper .dataTables_processing, 
      .dataTables_wrapper .dataTables_paginate { color: #e0e0e0 !important; }
      
      .dataTable thead th { background-color: #333 !important; color: #fff !important; border-color: #555 !important; }
      .dataTable tbody tr { background-color: #2c2c2c !important; color: #e0e0e0 !important; }
      .dataTable tbody tr:hover { background-color: #383838 !important; }
      .dataTable tbody tr.selected { background-color: #375a7f !important; color: #fff !important; }
      
      /* Inputs et Selects */
      input, select, .selectize-input { background-color: #333 !important; color: #fff !important; border-color: #555 !important; }
      .selectize-dropdown { background-color: #333 !important; color: #fff !important; }
      .selectize-dropdown-content .option { color: #e0e0e0 !important; }
      .selectize-dropdown-content .option.highlight { background-color: #375a7f !important; }
      
      /* Scrollbars (Webkit) */
      ::-webkit-scrollbar { width: 10px; }
      ::-webkit-scrollbar-track { background: #1a1a1a; }
      ::-webkit-scrollbar-thumb { background: #555; border-radius: 5px; }
      ::-webkit-scrollbar-thumb:hover { background: #777; }
    "))
  )
)


server <- function(input, output, session) {
  
  # ==============================================================================
  # 0. TEXTES STATIQUES
  # ==============================================================================
  output$description <- renderText("This application analyzes BirdNET detection results. Select your folders on the left, then explore detections, spectrograms, and phenology graphs.")
  output$contributions <- renderText("Contributions: Amandine Serrurier, Jean-Nicolas Pradervand & Christophe Sahli\nSwiss Ornithological Institute")
  output$license <- renderText("MIT License © 2026 vogelwarte.ch")
  

  # ==============================================================================
  # 1. GESTION DES DOSSIERS (ShinyFiles)
  # ==============================================================================
  shinyDirChoose(
    input, "dir1", 
    roots = c(home = roots_home), 
    session = session,
    # restrictions = c(home = roots_home)
  )
  
  shinyDirChoose(
    input, "dir2", 
    roots = c(home = roots_home), 
    session = session,
    # restrictions = c(home = roots_home)
  )
  
  # Variables chemins
  dir1 <- reactiveVal(NULL)
  dir2 <- reactiveVal(NULL)

  parse_safe_path <- function(input_val, roots) {
    if (is.null(input_val) || !is.list(input_val)) return(NULL)
    tryCatch({
      path <- parseDirPath(roots, input_val)
      if (is.character(path) && length(path) > 0 && path != "" && file.exists(path)) {
        return(path)
      }
    }, error = function(e) {
      showNotification(paste("Erreur chemin:", e$message), type = "warning")
    })
    return(NULL)
  }
  

  dir1_path <- reactive({
    req(input$dir1)
    parse_safe_path(input$dir1, c(home = roots_home))
  })
  
  dir2_path <- reactive({
    req(input$dir2)
    parse_safe_path(input$dir2, c(home = roots_home))
  })
  

  # Affichage chemins
  output$dir1_path <- renderText({
    p <- dir1_path()
    if (!is.null(p)) paste("Results:", basename(p)) else "No folder selected"
  })
  
  output$dir2_path <- renderText({
    p <- dir2_path()
    if (!is.null(p)) paste("Audio:", basename(p)) else "No folder selected"
  })
  
  # ==============================================================================
  # 2. CHARGEMENT DES DONNÉES
  # ==============================================================================

  DT_reac <- reactiveVal(NULL)
  Sp_List <- reactiveVal(NULL)
  complete_audio_path <- reactiveVal(NULL)
  Date_DT_min <- reactiveVal(NULL)
  Date_DT_max <- reactiveVal(NULL)
  
  observeEvent(input$start, {
    req(input$dir1)
    
    # 1. Initialisation : On vide les anciennes données et on lance la barre
    DT_reac(NULL)
    Sp_List(NULL)
    
    # Utilisation de withProgress pour créer un contexte de chargement
    # message: Titre de la boite de dialogue
    # detail: Texte secondaire (optionnel)
    withProgress(message = "Loading BirdNET Data...",
                 detail = "Initializing file scan...",
                 value = 0,
                 min = 0, max = 1, {
                   
                   # Étape 1 : Vérification fonction (rapide) -> +10%
                   incProgress(0.1, detail = "Checking functions...")
                   
                   if (!exists("load_selection_tables")) {
                     showNotification("Error: Function 'load_selection_tables' not found.", type = "error")
                     return()
                   }
                   
                   # Étape 2 : Lancement du chargement lourd -> +40% pendant le processus
                   incProgress(0.1, detail = "Scanning folders and reading tables...")
                   
                   tryCatch({
                     # La fonction load_selection_tables s'exécute ici
                     # Si cette fonction est elle-même longue, la barre restera à 20% jusqu'à la fin de cette ligne
                     DT <- load_selection_tables(
                       dir1 = dir1_path(), 
                       dir2 = dir2_path(), 
                       compiled = input$Compiled_F, 
                       device_tz = input$device_tz,      # Ex: "America/Anchorage" (choisi par l'utilisateur)
                       deployment_tz = input$deployment_tz # Ex: "Europe/Zurich" (choisi par l'utilisateur)
                     )
                     
                     # Étape 3 : Données lues, traitement en cours -> +70%
                     incProgress(0.4, detail = "Processing dates and species lists...")
                     
                     if (is.null(DT) || nrow(DT) == 0) {
                       showNotification("No data found or empty tables.", type = "error")
                       output$status <- renderText("No data loaded.")
                       return()
                     }
                     
                     # Étape 4 : Mise à jour des réactifs -> +90%
                     incProgress(0.2, detail = "Updating interface...")
                     
                     complete_audio_path(DT$Begin_Path)
                     Date_DT_min(min(DT$Date, na.rm = TRUE))
                     Date_DT_max(max(DT$Date, na.rm = TRUE))
                     
                     updateDateRangeInput(session, "filter_dates", 
                                          start = Date_DT_min(), end = Date_DT_max(),
                                          min = Date_DT_min(), max = Date_DT_max())
                     
                     Sp_List_DT <- sort(unique(DT$Common.Name))
                     Sp_List(Sp_List_DT)
                     
                     species_list_pheno <- c(Sp_List_DT, "All species")
                     
                     if (!is.null(input$species_pheno)) {
                       updateSelectizeInput(session, "species_pheno", choices = species_list_pheno, 
                                            selected = "All species", server = TRUE)
                     }
                     
                     # Sauvegarde finale
                     DT_reac(DT)
                     
                     # Étape 5 : Terminé -> 100% et fermeture automatique
                     incProgress(0.1, detail = "Done!")
                     output$status <- renderText(paste("Loaded:", nrow(DT), "detections |", length(Sp_List_DT), "species"))
                     
                   }, error = function(e) {
                     # En cas d'erreur, on ferme la barre et on affiche l'erreur
                     showNotification(paste("Error loading data:", e$message), type = "error")
                     output$status <- renderText("Error loading data.")
                     # La barre se fermera toute seule à la sortie du withProgress
                   })
                 }) # Fin du withProgress
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

  # ==============================================================================
  # 3. FILTRAGE RÉACTIF
  # ==============================================================================
  DT_filtered <- reactiveVal(NULL)
  last_species_list <- reactiveVal(NULL)
  
  filtered_confid <- reactive(input$Confid) %>% debounce(300)
  filtered_dates <- reactive(input$filter_dates) %>% debounce(300)
  
  observe({
    tryCatch({
      
      # 1. Vérification des prérequis
      req(DT_reac())
      DT_source <- DT_reac()
      n_rows <- nrow(DT_source)
      
      # Si le tableau source est vide, on arrête ici proprement
      if (n_rows == 0) {
        DT_filtered(DT_source)
        return()
      }
      
      # Récupération des inputs
      conf_val <- input$Confid
      dates_val <- input$filter_dates
      hour_range <- input$Hour_Range
      species_val <- input$species
      
      
      # Vérification de validité des inputs (doivent être des vecteurs de 2)
      # Si pas valides, on retourne TOUT le dataset sans filtrer
      if (is.null(conf_val) || length(conf_val) != 2 || 
          is.null(dates_val) || length(dates_val) != 2 || 
          is.null(hour_range) || length(hour_range) != 2) {
        DT_filtered(DT_source) # Retourne tout
        return()
      }
      
      DT_source <- DT_source[recorder %in% input$selected_recorders]
      
      
      # 2. Initialisation des vecteurs logiques (TRUE = on garde la ligne)
      condition_conf <- rep(TRUE, n_rows)
      condition_date <- rep(TRUE, n_rows)
      condition_hour <- rep(TRUE, n_rows)
      
      alerte_date <- FALSE
      alerte_hour <- FALSE
      message_alerte <- ""
      
      # --- Condition Confiance ---
      # On suppose que la colonne Confidence existe toujours si DT_source est valide
      if ("Confidence" %in% names(DT_source)) {
        condition_conf <- DT_source$Confidence >= conf_val[1] & 
          DT_source$Confidence <= conf_val[2]
      }
      
      # --- Condition Date ---
      if ("Date" %in% names(DT_source)) {
        dates_valides <- !is.na(DT_source$Date)
        
        if (sum(dates_valides) == 0) {
          # CAS : Aucune date valide -> On ne peut pas filtrer
          alerte_date <- TRUE
          condition_date <- rep(TRUE, n_rows) # On garde TOUT (fallback demandé)
        } else {
          # CAS : Dates valides présentes -> On applique le filtre
          condition_date <- dates_valides & 
            DT_source$Date >= dates_val[1] & 
            DT_source$Date <= dates_val[2]
        }
      } else {
        # Colonne absente
        alerte_date <- TRUE 
      }
      
      # --- Condition Hour ---
      if ("Hour" %in% names(DT_source)) {
        hours_valides <- !is.na(DT_source$Hour)
        
        if (sum(hours_valides) == 0) {
          # CAS : Aucune heure valide -> On ne peut pas filtrer
          alerte_hour <- TRUE
          condition_hour <- rep(TRUE, n_rows) # On garde TOUT
        } else {
          # CAS : Heures valides présentes -> On applique le filtre
          condition_hour <- hours_valides & 
            DT_source$Hour >= hour_range[1] & 
            DT_source$Hour <= hour_range[2]
        }
      } else {
        # Colonne absente
        alerte_hour <- TRUE
      }
      
      # 3. Gestion des alertes (Méthode robuste sans conflit de nom)
      if (alerte_date || alerte_hour) {
        if (alerte_date) message_alerte <- "Date invalide/absente. "
        if (alerte_hour) message_alerte <- paste0(message_alerte, "Heure invalide/absente. ")
        
        msg_final <- paste0("Filtrage partiel : ", message_alerte, "-> Affichage de toutes les données disponibles.")
        
        # Affichage dans la console R (toujours visible)
        cat("INFO FILTRE:", msg_final, "\n")
        
        # Notification Shiny (méthode directe pour éviter conflits)
        showNotification(
          msg_final,
          type = "warning",
          duration = 5,
          id = "filter_warning"
        )
      } else {
        removeNotification("filter_warning", session = session)
      }
      
      # 4. Application du filtre global
      global_filter <- condition_conf & condition_date & condition_hour
      
      # Filtrage
      DT_filter <- DT_source[global_filter, ]
      
      # SÉCURITÉ SUPPLÉMENTAIRE : Si le filtre vide tout le tableau (0 ligne)
      # Et que ce n'était pas voulu (ex: intervalle trop strict), on pourrait prévenir.
      # Mais ici on respecte le filtre logique.
      
      # 5. Mise à jour de la liste des espèces
      if (nrow(DT_filter) > 0) {
        available_species <- sort(unique(DT_filter$Common.Name))
      } else {
        available_species <- character(0)
      }
      
      current_selectized_list <- last_species_list()
      
      if (is.null(current_selectized_list) || !identical(current_selectized_list, available_species)) {
        updateSelectizeInput(
          session, 
          "species", 
          choices = available_species, 
          selected = if (!is.null(species_val) && species_val %in% available_species) species_val else NULL,
          server = TRUE
        )
        last_species_list(available_species)
      }
      
      # 6. Filtre Espèce
      if (!is.null(species_val) && species_val != "" && species_val %in% available_species) {
        DT_filter <- DT_filter[Common.Name == species_val, ]
      }
      
      # 7. Sauvegarde (Même si vide, on sauvegarde l'objet)
      DT_filtered(DT_filter)
      
    }, error = function(e) {
      # En cas d'erreur CRITIQUE (bug de code), on retourne TOUT le dataset pour ne pas bloquer l'app
      cat("ERREUR CRITIQUE DANS LE FILTRAGE:", e$message, "\n")
      showNotification(paste("Erreur technique (données complètes affichées):", e$message), type = "error", duration = 10)
      
      if (exists("DT_source")) {
        DT_filtered(DT_source) # Fallback : on montre tout
      } else {
        DT_filtered(NULL)
      }
    })
  })
  # ==============================================================================
  # 4. VARIABLES RÉACTIVES AUDIO (Déclarations uniques)
  # ==============================================================================
  selected_file <- reactiveVal(NULL)
  Name_reac <- reactiveVal(NULL)
  Date_reac <- reactiveVal(NULL)
  val_begin_start <- reactiveVal(NULL)
  val_begin_end <- reactiveVal(NULL)
  audio_error_msg <- reactiveVal(NULL)

  
  segment_audio <- reactiveVal(NULL)
  Samp_rate_stock <- reactiveVal(NULL)  
  wav_segment_stock_export <- reactiveVal(NULL) 
  bit_depth_reac <- reactiveVal(NULL) 
  audio_file_path <- reactiveVal(NULL)
  
  # ==============================================================================
  # 5. TRAITEMENT AUDIO PRINCIPAL (Déclenché par Sélection OU Highpass)
  # ==============================================================================
 
  # --- OBSERVATEUR NETTOYEUR (CRUCIAL) ---
  # Dès que le slider bouge, on EFFACE les données audio actuelles.
  # Cela force le graphique à se vider et empêche l'affichage de l'ancien résultat.
  observeEvent(input$Highpass, {
    # On met à NULL pour signaler "en cours de recalcul"
    segment_audio(NULL)
    audio_error_msg(NULL) 
    # Optionnel : montrer un message "Filtrage..."
    # showNotification("Filtrage en cours...", duration = 1)
  }, ignoreNULL = FALSE, ignoreInit = FALSE)
  
  
  
  observeEvent(list(input$summary_rows_selected, input$Highpass), {
    # 1. Prérequis stricts
    req(DT_filtered())
    req(input$summary_rows_selected)
    
    # 2. Récupération des données de la ligne
    selected_row <- input$summary_rows_selected[1]
    dataset_chosen <- DT_filtered()
    # dataset_chosen <- DT # debug
    # selected_row <- 1
    if (selected_row > nrow(dataset_chosen) || selected_row < 1) return()
    
    # 3. Mise à jour des métadonnées (Texte, Dates, etc.)
    # On ne fait ça QUE si la ligne change (optimisation), mais c'est inoffensif de le refaire
    file_path_val <- dataset_chosen$True_Location[selected_row]
    selected_file(file_path_val)
    Name_reac(dataset_chosen$Common.Name[selected_row])
    

    if (input$Compiled_F) {
      val_begin_start(dataset_chosen$File.Offset..s.[selected_row])
      val_begin_end(dataset_chosen$File.Offset..s.[selected_row] + 3)
    } else {
      val_begin_start(dataset_chosen$Begin.Time..s.[selected_row])
      val_begin_end(dataset_chosen$End.Time..s.[selected_row])
    }
    
    # Date
    Date_temp <- dataset_chosen$Date[selected_row]
    if (is.na(Date_temp)) Date_temp <- dataset_chosen$file_name_basic[selected_row]
    else Date_temp <- paste0(Date_temp, "_", dataset_chosen$Hour[selected_row], "h", 
                             format(dataset_chosen$Min[selected_row], fmt="%02d"), "min")
    Date_reac(Date_temp)
    
    output$selected_audio <- renderText({
      if(is.null(file_path_val) || !nzchar(file_path_val)) "No file" else paste("File:", basename(file_path_val))
    })

    
    # 4. VÉRIFICATION FICHIER
    if (is.null(file_path_val) || !file.exists(file_path_val)) {
      audio_error_msg("Fichier introuvable")
      segment_audio(NULL)
      return()
    }
    
    # 5. LECTURE DIRECTE ET IMMÉDIATE DU HIGHPASS
    # C'est l'étape CRUCIALE. On lit input$Highpass ICI, dans le contexte actif.
    # On ne passe par aucune variable intermédiaire.
    current_hp_val <- input$Highpass
    
    # Debug : Décommentez la ligne ci-dessous pour voir dans la console R quelle valeur est lue
    # cat(">>> TRAITEMENT AVEC HIGHPASS =", current_hp_val, "\n")
    
    if (!exists("prepare_for_spectro")) {
      audio_error_msg("Fonction manquante")
      segment_audio(NULL)
      return()
    }
    
    # 6. EXÉCUTION DU TRAITEMENT
    tryCatch({
      # On passe current_hp_val directement
      list_of_res <- prepare_for_spectro(
        file_path = file_path_val, 
        Begin_Start = val_begin_start(), 
        Begin_End = val_begin_end(), 
        highpass = current_hp_val 
      )
      
      # 7. MISE À JOUR ATOMIQUE
      # Si on arrive ici, c'est que le traitement a réussi avec la NOUVELLE valeur
      segment_audio(list_of_res[[1]])
      Samp_rate_stock(list_of_res[[2]])
      wav_segment_stock_export(list_of_res[[3]])
      bit_depth_reac(list_of_res[[4]])
      audio_error_msg(NULL)
      
    }, error = function(e) {
      audio_error_msg(paste("Erreur:", e$message))
      segment_audio(NULL)
      # Affiche l'erreur dans la console pour debug
      # cat("ERREUR TRAITEMENT:", e$message, "\n")
    })
  }, ignoreNULL = FALSE, ignoreInit = FALSE)
  

  # ==============================================================================
  # 6. GÉNÉRATION FICHIER TEMPORAIRE (Lecteur Audio)
  # ==============================================================================
  observeEvent(segment_audio(), {
    req(segment_audio())
    req(Samp_rate_stock())
    
    old_file <- audio_file_path()
    if (!is.null(old_file) && file.exists(old_file)) {
      tryCatch(file.remove(old_file), error = function(e) invisible(NULL))
    }
    
    wav_segment <- segment_audio()
    if (!inherits(wav_segment, "Wave")) {
      wav_segment <- Wave(wav_segment, samp.rate = Samp_rate_stock(), bit = bit_depth_reac())
    }
    
    wav_mono <- tryCatch(mono(wav_segment, which = "left"), error = function(e) wav_segment)
    
    if (!is.integer(wav_mono@left)) {
      max_val <- max(abs(wav_mono@left))
      if (max_val > 0) {
        wav_mono@left <- round((wav_mono@left / max_val) * 32767)
        wav_mono@bit <- 16
      }
    }
    
    new_temp_file <- tempfile(fileext = ".wav")
    success <- tryCatch({ writeWave(wav_mono, new_temp_file); TRUE }, error = function(e) FALSE)
    
    if (success && file.exists(new_temp_file)) {
      audio_file_path(new_temp_file)
    }
  })
  
  # ==============================================================================
  # 7. AFFICHAGE LECTEUR & NETTOYAGE
  # ==============================================================================
  output$audio_player <- renderUI({
    file_path <- audio_file_path()
    if (is.null(file_path) || !file.exists(file_path)) {
      return(tags$p("En attente de l'audio...", style = "color: #888;"))
    }
    
    wav_base64 <- tryCatch(base64enc::base64encode(file_path), error = function(e) NULL)
    if (is.null(wav_base64)) return(tags$p("Erreur encodage", style="color:red;"))
    
    tags$audio(src = paste0("data:audio/wav;base64,", wav_base64), controls = TRUE, style = "width: 100%;")
  })
  
  session$onSessionEnded(function() {
    f <- audio_file_path()
    if (!is.null(f) && file.exists(f)) file.remove(f)
  })
  
  # ==============================================================================
  # 8. OUVERTURE FICHIER EXTERNE
  # ==============================================================================
  observeEvent(input$open_audio, {
    req(selected_file())
    file_path <- selected_file()
    
    if (!file.exists(file_path)) {
      showModal(modalDialog("Error: Audio file not found on disk.", easyClose = TRUE))
      return()
    }
    
    os_type <- Sys.info()["sysname"]
    cmd <- NULL
    
    if (os_type == "Windows") {
      cmd <- paste0('cmd.exe /c start "" ', shQuote(file_path))
    } else if (os_type == "Darwin") {
      cmd <- paste0('open ', shQuote(file_path))
    } else if (os_type == "Linux") {
      if (system2("which", args = "xdg-open", stdout = TRUE, stderr = TRUE) == "") {
        showModal(modalDialog("Command 'xdg-open' not found.", easyClose = TRUE))
        return()
      }
      cmd <- paste0('xdg-open ', shQuote(file_path))
    } else {
      showModal(modalDialog(paste("Unsupported OS:", os_type), easyClose = TRUE))
      return()
    }
    
    exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (exit_code != 0) {
      showNotification(paste("Failed to open file. Exit code:", exit_code), type = "warning")
    }
  })
  
  # ==============================================================================
  # 9. TABLEAU DE DONNÉES
  # ==============================================================================
  output$summary <- renderDT({

    req(DT_filtered())
    dataset_chosen <- DT_filtered()
    
    if (nrow(dataset_chosen) == 0) {
      return(datatable(data.frame(Message = "Aucune donnée ne correspond aux filtres"), 
                       options = list(dom = 't')))
    }
    
    cols_to_show <- intersect(c("Common.Name", "Date", "Hour", "Min", "DateTime_Display",
                                "Begin.Time..s.", "End.Time..s.", "Confidence"), 
                              names(dataset_chosen))
    
    data_to_show <- dataset_chosen[, ..cols_to_show]
    

    if ("Begin.Time..s." %in% names(data_to_show)) 
      data_to_show[, Begin.Time..s. := round(Begin.Time..s., 0)]
    if ("End.Time..s." %in% names(data_to_show)) 
      data_to_show[, End.Time..s. := round(End.Time..s., 0)]
    
    datatable(data_to_show, 
              escape = FALSE, 
              selection = 'single', 
              options = list(pageLength = 15, scrollX = TRUE, autoWidth = TRUE))
  })
  
  # ==============================================================================
  # 10. EXPORTS
  # ==============================================================================
  output$download <- downloadHandler(
    filename = function() { paste0(Sys.Date(), "_BirdNet_Data.csv") },
    content = function(fname) { 
      req(DT_filtered())
      write.csv(DT_filtered(), fname, row.names = FALSE) 
    }
  )
  

  output$download_audio <- downloadHandler(
    filename = function() { paste0(clean_text(Name_reac()), "_segment_", Date_reac(), ".wav") },
    content = function(file) {
      req(wav_segment_stock_export())
      if (exists("export_selected_audio")) {
        wav_out <- export_selected_audio(wav = wav_segment_stock_export(), 
                                         Samp_rate = Samp_rate_stock(), 
                                         bit = bit_depth_reac(), 
                                         highpass = input$Highpass)
        writeWave(wav_out, file)
      } else {
        writeWave(wav_segment_stock_export(), file)
      }
    }
  )
  
  output$download_audio_mp3 <- downloadHandler(
    filename = function() { 

      paste0(clean_text(Name_reac()), "_segment_", Date_reac(), ".mp3") 
    },
    content = function(file) {
      req(wav_segment_stock_export())
      if (exists("export_selected_audio")) {
        wav_out <- export_selected_audio(wav = wav_segment_stock_export(), 
                                         Samp_rate = Samp_rate_stock(), 
                                         bit = bit_depth_reac(), 
                                         highpass = input$Highpass)
        writeWave(wav_out, file)
      } else {
        writeWave(wav_segment_stock_export(), file)
      }
    }
  )
  
  # ==============================================================================
  # 11. SPECTROGRAMME & ERREURS
  # ==============================================================================
  output$Spectro <- renderPlot({
    req(segment_audio())
    req(Samp_rate_stock())
    
    # --- LIGNE MAGIQUE POUR FORCER LA MISE À JOUR ---
    # On lit input$Highpass ici. Si ça change, renderPlot se relance OBLIGATOIREMENT.
    # On l'assigne à une variable inutile juste pour créer la dépendance.
    force(input$Highpass) 
    # -----------------------------------------------
    
    
    tryCatch({
      spectro(
        segment_audio(), 
        norm = TRUE, 
        flim = (input$Freq_Spectro/1000), 
        f = Samp_rate_stock(),
        wl = 1024,
        ovlp = 50,
        #dB = TRUE,             # <--- AJOUT : Affichage en Décibels (meilleur contraste)
        collevels = seq(-80, 0, 1), # <--- AJOUT : Échelle de -80dB (bruit) à 0dB (max)
        palette = reverse.gray.colors.1, 
        noisereduction = 1, 
        main = paste0(Name_reac(), " | ", Samp_rate_stock(), " Hz"),
        oma = c(1, 1, 3, 1)
      )
    }, error = function(e) {
      plot(1, type="n", axes=FALSE, main="Erreur Spectro")
      text(0.5, 0.5, e$message, col="red")
    })
  })
  
  output$Spectro_error <- renderUI({
    err <- audio_error_msg()
    if (!is.null(err) && err != "") {
      tags$div(style = "color: #ff6b6b; text-align: center; padding: 20px;", 
               h4("Erreur de chargement"), p(err))
    }
  })
  
} # Fin du server

# --- 6. LANCEMENT ---
shinyApp(ui = ui, server = server)



