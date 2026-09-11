pheno_matrix <- function(Voc, SP = "All species", Unit, Confidence1, sunrise, LAT, LONG, TimeZone, xlim_plot = NULL, nocturnal = FALSE) {
  
  # Chargement explicite pour éviter les erreurs de namespace
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Package data.table required")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package ggplot2 required")
  if (!requireNamespace("suncalc", quietly = TRUE)) stop("Package suncalc required")
  if (!requireNamespace("scales", quietly = TRUE)) stop("Package scales required")
  if (!requireNamespace("lubridate", quietly = TRUE)) stop("Package lubridate required")

  # --- 1. OPTIMISATION DU FILTRAGE ---
  # Conversion sécurisée en data.table
  if (!is.data.table(Voc)) {
    dt <- as.data.table(Voc)
  } else {
    dt <- Voc
  }
  
  # # FILTRAGE PRIORITAIRE PAR ESPÈCE
  # if (!"All species" %in% SP) {
  #   dt <- dt[Common.Name %in% SP]
  # }
  
  # FILTRAGE PAR CONFIANCE
  #dt <- dt[Confidence >= Confidence1]
  
  if (nrow(dt) == 0) {
    return("Error: species not present in the dataset or confidence threshold set too high for the species selected")
  }
  
  # Titre
  if (length(SP) > 1) {
    SP_title <- paste0("Multiple detected species considered, n = ", length(unique(dt$Common.Name)), " species")
  } else {
    SP_title <- if (SP == "All species") {
      paste0("All detected species considered, n = ", length(unique(dt$Common.Name)), " species")
    } else {
      SP
    }
  }
  
  # --- 2. PRÉPARATION TEMPS (CORRIGÉ FUSEAU HORAIRE) ---
  SeqMin <- seq(0, (1440 - Unit), by = Unit)
  
  # Calcul temps en secondes depuis minuit LOCAL (et non UTC)
  if (inherits(dt$Start_segment, "POSIXct")) {
    # Extraction explicite des composants H, M, S dans le fuseau de l'objet (deployment_tz)
    # format() respecte l'attribut 'tzone', contrairement à as.numeric()
    dt[, `:=`(
      h_loc = as.numeric(format(Start_segment, "%H")),
      m_loc = as.numeric(format(Start_segment, "%M")),
      s_loc = as.numeric(format(Start_segment, "%S"))
    )]
    
    # Calcul vectoriel : (Heures * 3600) + (Minutes * 60) + Secondes
    dt[, time_numeric := (h_loc * 3600) + (m_loc * 60) + s_loc]
    
    # Nettoyage des colonnes temporaires pour garder le tableau propre
    dt[, c("h_loc", "m_loc", "s_loc") := NULL]
    
  } else {
    # Fallback sécurisé (si le format n'est pas POSIXct)
    # Note : as_hms peut avoir ses propres comportements de timezone, à utiliser en dernier recours
    dt[, time_numeric := as.numeric(as_hms(Start_segment))]
  }
  
  # Bornes pour le découpage en intervalles
  SeqMin_numeric <- SeqMin * 60
  max_time <- max(dt$time_numeric, na.rm = TRUE)
  breaks <- unique(c(SeqMin_numeric, max_time + 1))
  breaks <- sort(breaks) # Important pour la fonction cut()
  
  # Dates et Année (Basés sur l'heure LOCALE)
  # On extrait l'année du premier élément pour construire le calendrier de référence
  year1 <- as.numeric(format(dt$Start_segment[1], "%Y"))
  
  # Calcul du Day Of Year (1 à 366) dans le fuseau local
  dt[, DOY := as.integer(format(Start_segment, "%j"))]
  
  # Création du calendrier de référence pour la matrice
  SeqDate <- data.frame(
    date = seq.Date(as.Date(paste0(year1, "-01-01")), as.Date(paste0(year1, "-12-31")), by = "day"),
    stringsAsFactors = FALSE
  )
  SeqDate$DOY <- as.integer(format(SeqDate$date, "%j"))
  
  
  
  
  # --- 3. MATRICE D'AGRÉGATION ---
  dt[, interval := cut(time_numeric, breaks = breaks, right = FALSE, include.lowest = TRUE, labels = FALSE)]
  
  # Agrégation
  counts <- dt[, .N, by = .(interval, DOY)]
  counts <- counts[!is.na(interval)] # Nettoyage
  
  # Initialisation Matrice
  n_time <- length(SeqMin)
  n_days <- nrow(SeqDate)
  VocMatrix <- matrix(0L, nrow = n_time, ncol = n_days, dimnames = list(NULL, SeqDate$DOY))
  
  # Remplissage Vectorisé (Sécurisé)
  valid_rows <- counts$interval[counts$interval <= n_time]
  valid_cols_idx <- match(counts$DOY[counts$interval <= n_time], SeqDate$DOY)
  valid_vals <- counts$N[counts$interval <= n_time]
  
  if (length(valid_rows) > 0 && !any(is.na(valid_cols_idx))) {
    VocMatrix[cbind(valid_rows, valid_cols_idx)] <- valid_vals
  }
  
  # --- 4. FORMAT LONG POUR GGPLOT ---
  matrix_df <- as.data.frame(VocMatrix)
  matrix_df$seqY <- SeqMin / 60
  
  # Conversion EXPLICITE en data.table AVANT le melt
  matrix_dt <- as.data.table(matrix_df)
  
  # Melt natif data.table
  matrix_long <- melt(matrix_dt, id.vars = "seqY", variable.name = "DOY", value.name = "Vocs")
  
  # Nettoyage et calculs
  matrix_long[, DOY := as.integer(as.character(DOY))]
  matrix_long[, date_vocmatrix := as.Date(DOY - 1, origin = paste0(year1, "-01-01"))]
  matrix_long[, time := format(as.POSIXct("00:00:00", format="%H:%M:%S") + (seqY * 3600), format="%H:%M")]
  
  # Filtre xlim_plot
  if (!is.null(xlim_plot)) {
    SeqDateDOY <- SeqDate[SeqDate$date >= xlim_plot[1] & SeqDate$date <= xlim_plot[2], ]$DOY
    matrix_long <- matrix_long[DOY %in% SeqDateDOY]
  }
  
  # Conversion pour ggplot
  ref_date <- as.numeric(as.Date("2000-01-01")) * 86400
  matrix_long[, time_ct := as.POSIXct(time, format="%H:%M") + ref_date]
  
  # --- 5. CONSTRUCTION DU PLOT ---
  
  # Vérification de sécurité
  if (nrow(matrix_long) == 0) {
    return("Error: No data to plot after processing.")
  }
  
  if (!nocturnal) {
    # === MODE DIURNE (AXE 0-24H SIMPLE) ===
    
    # 1. Calcul de l'heure décimale (0 à 24)
    # On part de la colonne 'time' originale (format HH:MM)
    matrix_long[, hour_base := as.numeric(format(as.POSIXct(time, format="%H:%M"), format="%H")) + 
                  (as.numeric(format(as.POSIXct(time, format="%H:%M"), format="%M")) / 60)]
    
    if (all(is.na(matrix_long$hour_base))) {
      stop("Error: All time values are NA in diurnal mode.")
    }
    
    # Pour le diurne, hour_num = hour_base (pas de décalage)
    matrix_long[, hour_num := hour_base]
    
    # Nettoyage
    matrix_long <- matrix_long[!is.na(hour_num) & !is.na(date_vocmatrix) & !is.na(Vocs)]
    
    # 2. Données Solaires (Conversion en heures décimales 0-24)
    SunDF <- NULL
    if (sunrise) {
      unique_dates <- unique(matrix_long$date_vocmatrix)
      if (length(unique_dates) > 0) {
        tryCatch({
          Sun <- getSunlightTimes(date = unique_dates, lat = LAT, lon = LONG, keep = c("sunrise", "sunset", "dawn", "dusk"), tz = TimeZone)
          
          # Helper simple 0-24
          to_hour_num <- function(posix_time) {
            if (is.na(posix_time)) return(NA)
            h <- as.numeric(format(posix_time, "%H"))
            m <- as.numeric(format(posix_time, "%M"))
            return(h + m/60)
          }
          
          SunDF <- data.frame(
            date_vocmatrix = unique_dates,
            ymin_dawn = sapply(Sun$dawn, to_hour_num),
            ymax_dusk = sapply(Sun$dusk, to_hour_num), # Note: dusk est souvent > dawn, mais pour geom_ribbon on veut ymin/ymax corrects
            ymin_sunrise = sapply(Sun$sunrise, to_hour_num),
            ymax_sunset = sapply(Sun$sunset, to_hour_num),
            stringsAsFactors = FALSE
          )
          
          # Correction pour geom_ribbon : ymin doit être < ymax
          # Pour l'aube (dawn -> sunrise) : ymin=dawn, ymax=sunrise
          # Pour le jour (sunrise -> sunset) : pas de ribbon (ou blanc)
          # Pour le crépuscule (sunset -> dusk) : ymin=sunset, ymax=dusk
          # On simplifie ici en gardant les rubans "nuit" et "crépuscule/aube"
          
          # Réorganisation pour correspondre à la logique de plot ci-dessous
          SunDF$ymin_morning <- pmin(SunDF$ymin_dawn, SunDF$ymin_sunrise)
          SunDF$ymax_morning <- pmax(SunDF$ymin_dawn, SunDF$ymin_sunrise)
          
          SunDF$ymin_evening <- pmin(SunDF$ymax_sunset, SunDF$ymax_dusk)
          SunDF$ymax_evening <- pmax(SunDF$ymax_sunset, SunDF$ymax_dusk)
          
          SunDF <- SunDF[complete.cases(SunDF), ]
          if (nrow(SunDF) == 0) SunDF <- NULL
          
        }, error = function(e) {
          warning("Sunlight calculation failed: ", e$message)
        })
      }
    }
    
    # 3. Construction du graphique
    a <- ggplot()
    
    # --- COUCHE 1 : SOLEIL (FOND) ---
    if (!is.null(SunDF) && nrow(SunDF) > 0 && sunrise) {
      a <- a +
        # Nuit (Orange) : 0 -> Dawn (si dawn > 0)
        geom_ribbon(data = SunDF, aes(x = date_vocmatrix, ymin = 0, ymax = ymin_morning), fill = "dodgerblue4", alpha = 0.4) +
        # Nuit (Orange) : Dusk -> 24
        geom_ribbon(data = SunDF, aes(x = date_vocmatrix, ymin = ymax_evening, ymax = 24), fill = "dodgerblue4", alpha = 0.4) +
        # Aube (Bleu) : Dawn -> Sunrise
        geom_ribbon(data = SunDF, aes(x = date_vocmatrix, ymin = ymin_morning, ymax = ymax_morning), fill = "#CD853F", alpha = 0.4) +
        # Crépuscule (Bleu) : Sunset -> Dusk
        geom_ribbon(data = SunDF, aes(x = date_vocmatrix, ymin = ymin_evening, ymax = ymax_evening), fill = "#CD853F", alpha = 0.4)
    }

    # --- COUCHE 2 : TUILES (DESSUS) ---
    matrix_plot <- matrix_long[!is.na(hour_num) & !is.na(date_vocmatrix) & Vocs > 0]
    
    if (nrow(matrix_plot) > 0) {
      df_plot <- as.data.frame(matrix_plot)
      a <- a +
        geom_tile(
          data = df_plot, 
          aes(x = date_vocmatrix, y = hour_num, fill = Vocs), 
          color = NA, alpha = 1
        ) +
        scale_fill_gradientn(
          colours = c("#FCFFA4FF", "#F98C0AFF", "#BB3754FF", "#56106EFF", "#000004FF"),
          na.value = "transparent",
          values = scales::rescale(c(0.01, 0.2, 0.4, 0.6, 1))
        )
    } else {
      warning("Diurnal Plot: No detections > 0 found.")
      a <- a + scale_fill_gradientn(
        colours = c("#FCFFA4FF", "#F98C0AFF", "#BB3754FF", "#56106EFF", "#000004FF"),
        na.value = "transparent",
        values = scales::rescale(c(0.01, 0.2, 0.4, 0.6, 1)),
        guide = "colorbar"
      )
    }
    
    # --- FINITION : AXE Y 0-24H ---
    a <- a +
      scale_y_continuous(
        "Time",
        limits = c(0, 24),           # 0h à 24h
        breaks = seq(0, 24, by = 2), # Un trait tous les 2h
        labels = function(x) sprintf("%02d:00", as.integer(x)),
        expand = c(0, 0)
      ) +
      theme(
        panel.background = element_rect(fill = "white"), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white"),
        text = element_text(size = 18)
      ) +
      labs(
        fill = paste0("Number of\ndetections/", Unit, " mins"), 
        y = "Time", 
        x = "Date", 
        title = SP_title
      )
    
  } else {

    
    # === MODE NOCTURNE (CORRIGÉ : HEURES DU MATIN INCLUSES) ===
    
    today <- Sys.Date()
    today2 <- today + 1
    
    # 1. Calcul de l'heure de base (0-24)
    matrix_long[, hour_base := as.numeric(format(as.POSIXct(time, format="%H:%M"), format="%H")) + 
                  (as.numeric(format(as.POSIXct(time, format="%H:%M"), format="%M")) / 60)]
    
    if (all(is.na(matrix_long$hour_base))) {
      warning("Nocturnal Plot Error: Impossible de lire les heures.")
      return(NULL)
    }
    
    # 2. Création de la date graphique (pour l'axe X)
    # Si heure >= 12 -> Jour J (today)
    # Si heure < 12  -> Jour J+1 (today2) pour l'affichage
    matrix_long[, graph_date_str := ifelse(hour_base >= 12, as.character(today), as.character(today2))]
    matrix_long[, test_graph := as.POSIXct(paste(graph_date_str, time), format="%Y-%m-%d %H:%M", tz = "UTC")]
    
    # 3. Ajustement de l'axe X (date_vocmatrix)
    # On décale d'un jour vers l'arrière les points du matin (< 12h) pour les coller à la veille
    matrix_long[hour_base < 12, date_vocmatrix := date_vocmatrix - 1]
    
    # 4. CRÉATION DE L'AXE Y NUMÉRIQUE CONTINU (12h -> 36h)
    # C'est ICI que se fait la magie :
    # - Si heure >= 12 : on garde la valeur (ex: 13.0)
    # - Si heure < 12  : on ajoute 24 (ex: 02.0 devient 26.0)
    matrix_long[, hour_num := ifelse(hour_base >= 12, hour_base, hour_base + 24)]
    
    # 5. NETTOYAGE
    matrix_long <- matrix_long[!is.na(hour_num) & !is.na(date_vocmatrix) & !is.na(Vocs)]
    
    if (nrow(matrix_long) == 0) {
      warning("Nocturnal Plot Error: Aucune donnée valide.")
      return(NULL)
    }
    
    # 6. Données Solaires (Conversion en heures continues 12-36h)
    SunDF <- NULL
    if (sunrise) {
      unique_dates <- unique(matrix_long$date_vocmatrix)
      if (length(unique_dates) > 0) {
        tryCatch({
          Sun <- getSunlightTimes(date = unique_dates, lat = LAT, lon = LONG, keep = c("sunrise", "sunset", "dawn", "dusk"), tz = TimeZone)
          
          # Helper pour convertir en heure continue (12-36h)
          # Attention : sunrise/sunset sont sur des jours différents dans la logique graphique
          to_continuous_hour <- function(posix_time, is_next_day = FALSE) {
            if (is.na(posix_time)) return(NA)
            h <- as.numeric(format(posix_time, "%H"))
            m <- as.numeric(format(posix_time, "%M"))
            val <- h + m/60
            if (is_next_day) return(val + 24) # Si c'est le lendemain graphique, on ajoute 24
            return(val)
          }
          
          # Logique : 
          # Sunset/Dusk (Jour J) -> Reste 12-24
          # Sunrise/Dawn (Jour J+1) -> Devient 24-36
          
          SunDF <- data.frame(
            date_vocmatrix = unique_dates,
            ymin_sunset = sapply(as.POSIXct(paste(format(today, "%Y-%m-%d"), format(Sun$sunset, "%H:%M:%S")), format="%Y-%m-%d %H:%M:%S", tz = TimeZone), to_continuous_hour, is_next_day=FALSE),
            ymax_midnight = 24.0, 
            ymax_sunrise = sapply(as.POSIXct(paste(format(today2, "%Y-%m-%d"), format(Sun$sunrise, "%H:%M:%S")), format="%Y-%m-%d %H:%M:%S", tz = TimeZone), to_continuous_hour, is_next_day=TRUE),
            ymin_dusk = sapply(as.POSIXct(paste(format(today, "%Y-%m-%d"), format(Sun$dusk, "%H:%M:%S")), format="%Y-%m-%d %H:%M:%S", tz = TimeZone), to_continuous_hour, is_next_day=FALSE),
            ymax_dawn = sapply(as.POSIXct(paste(format(today2, "%Y-%m-%d"), format(Sun$dawn, "%H:%M:%S")), format="%Y-%m-%d %H:%M:%S", tz = TimeZone), to_continuous_hour, is_next_day=TRUE),
            stringsAsFactors = FALSE
          )
          
          SunDF <- SunDF[complete.cases(SunDF), ]
          if (nrow(SunDF) == 0) SunDF <- NULL
          
        }, error = function(e) {
          warning("Calcul solaire échoué : ", e$message)
        })
      }
    }
    
    # 7. Construction du graphique
    a <- ggplot()
    
    # --- COUCHE 1 : SOLEIL ---
    if (!is.null(SunDF) && nrow(SunDF) > 0 && sunrise) {
      a <- a +
        # Nuit (Orange) : Sunset (18h) -> Minuit (24h)
        geom_ribbon(data = SunDF, aes(x = date_vocmatrix, ymin = ymin_sunset, ymax = ymax_midnight), fill = "#CD853F", alpha = 0.4) +
        # Nuit (Orange) : Minuit (24h) -> Sunrise (ex: 30h = 06h)
        geom_ribbon(data = SunDF, aes(x = date_vocmatrix, ymin = 24, ymax = ymax_sunrise), fill = "#CD853F", alpha = 0.4) +
      # Crépuscule (Bleu) : Dusk -> Minuit
      geom_ribbon(data = SunDF, aes(x = date_vocmatrix, ymin = ymin_dusk, ymax = ymax_midnight), fill = "dodgerblue4", alpha = 0.4) +
        # Aube (Bleu) : Minuit -> Dawn
        geom_ribbon(data = SunDF, aes(x = date_vocmatrix, ymin = 24, ymax = ymax_dawn), fill = "dodgerblue4", alpha = 0.4)
    }
    
    # --- COUCHE 2 : TUILES ---
    matrix_plot <- matrix_long[!is.na(hour_num) & !is.na(date_vocmatrix) & Vocs > 0]
    
    if (nrow(matrix_plot) > 0) {
      df_plot <- as.data.frame(matrix_plot)
      a <- a +
        geom_tile(
          data = df_plot, 
          aes(x = date_vocmatrix, y = hour_num, fill = Vocs), 
          color = NA, alpha = 1
        ) +
        scale_fill_gradientn(
          colours = c("#FCFFA4FF", "#F98C0AFF", "#BB3754FF", "#56106EFF", "#000004FF"),
          na.value = "transparent",
          values = scales::rescale(c(0.01, 0.2, 0.4, 0.6, 1))
        )
    } else {
      warning("Nocturnal Plot: No detections > 0 found.")
      a <- a + scale_fill_gradientn(
        colours = c("#FCFFA4FF", "#F98C0AFF", "#BB3754FF", "#56106EFF", "#000004FF"),
        na.value = "transparent",
        values = scales::rescale(c(0.01, 0.2, 0.4, 0.6, 1)),
        guide = "colorbar"
      )
    }
    
    # --- FINITION : AXE Y PERSONNALISÉ (12h à 36h avec labels intelligents) ---
    a <- a +
      scale_y_continuous(
        "Time",
        limits = c(12, 36),          # 12h (midi) à 36h (06h du matin)
        breaks = seq(12, 36, by = 2), # Un trait tous les 2h : 12, 14, ..., 24, 26, ..., 36
        labels = function(x) {
          # Fonction personnalisée pour l'affichage :
          # Si x <= 24 -> affiche "12", "14", ..., "24"
          # Si x > 24  -> affiche "02", "04", ..., "12" (x - 24)
          real_hour <- ifelse(x >= 24, x - 24, x)
          # Formatage avec zéro devant si nécessaire
          sprintf("%02d:00", as.integer(real_hour))
        },
        expand = c(0, 0)
      ) +
      theme(
        panel.background = element_rect(fill = "white"), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white"),
        text = element_text(size = 18)
      ) +
      labs(
        fill = paste0("Number of\ndetections/", Unit, " mins"), 
        y = "Time", 
        x = "Date", 
        title = SP_title
      )
  }
  
  return(a)
}