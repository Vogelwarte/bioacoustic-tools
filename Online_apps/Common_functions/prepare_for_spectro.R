prepare_for_spectro <- function(file_path, Begin_Start, Begin_End, highpass) {
  
  # --- LIGNE DE DEBUG CRUCIALE ---
  message(">>> [INTERNE FONCTION] highpass reçu = ", highpass)
  # ------------------------------
  
  # --- 1. VÉRIFICATION ET PRÉPARATION ---
  if (is.null(file_path) || !file.exists(file_path)) {
    stop("Le fichier audio spécifié est introuvable.")
  }

  file_ext <- tolower(tools::file_ext(file_path))
  wav_file <- NULL
  
  # Paramètres de fenêtre
  margin <- 0.5 # Marge avant/après en secondes
  min_display_duration <- 3.0
  
  # --- 2. FONCTION UTILITAIRE POUR LIRE LA DURÉE SANS CHARGER LE FICHIER ---
  #path <- file_path
  get_audio_info <- function(path) {
    # Lecture de l'en-tête seulement
    info <- tryCatch({
      readWave(path, header = TRUE)
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(info)) return(NULL)
    
    # COMPATIBILITÉ TUNE R : 
    # Si info est une liste (nouvelles versions), on extrait les éléments nommés
    # Si info est un objet Wave (anciennes versions), on utilise @
    
    samp_rate <- NULL
    n_samples <- NULL
    bit_depth <- NULL
    
    if (inherits(info, "Wave")) {
      # Ancienne méthode (objet Wave)
      samp_rate <- info@samp.rate
      n_samples <- length(info@left) # Ou info@nl pour headeronly parfois
      bit_depth <- info@bit
    } else if (is.list(info)) {
      # Nouvelle méthode (liste) - Les noms peuvent varier selon la version
      # Souvent : samp.rate, bit, nsamp, etc.
      samp_rate <- info$sample.rate
      if (is.null(samp_rate)) samp_rate <- info[["sample.rate"]]
      
      bit_depth <- info$bits
      if (is.null(bit_depth)) bit_depth <- info[["bits"]]
      
      # Le nombre d'échantillons peut être dans 'nsamp' ou calculable
      n_samples <- info$samples
      if (is.null(n_samples)) n_samples <- info[["samples"]]
      
      # Si nsamp est manquant dans la liste headeronly, on doit parfois lire un peu plus
      # Mais normalement headeronly le fournit.
    }
    
    if (is.null(samp_rate) || is.null(n_samples)) return(NULL)
    
    # Calcul manuel de la durée (plus fiable que @duration qui n'existe pas toujours en headeronly)
    duration_sec <- n_samples / samp_rate
    
    return(list(samp.rate = samp_rate, bit = bit_depth, duration = duration_sec))
  }
  
  # --- 3. LOGIQUE DE LECTURE ---
  
  # Définition des temps cibles
  #Begin_Start = 0
  #Begin_End = 3
  req_start <- max(0, Begin_Start - margin)
  req_end <- Begin_End + margin
  
  file_info <- NULL
  
  if (file_ext %in% c("wav")) {
    # --- CAS WAV ---
    file_info <- get_audio_info(file_path)
    
    if (is.null(file_info)) {
      stop(paste(file_path, "Impossible de lire les informations du fichier WAV (corrompu ?)."))
    }
    
    Samp_rate <- file_info$samp.rate
    bit_depth_val <- file_info$bit
    full_duration <- file_info$duration
    
    # Calcul sécurisé des bornes de lecture
    read_to <- min(req_end, full_duration)
    read_from <- max(0, req_start)
    
    # Correction critique : Si from >= to
    if (read_from >= read_to) {
      read_to <- min(full_duration, read_from + min_display_duration)
      read_from <- max(0, read_to - min_display_duration)
    }
    
    # Lecture du segment
    tryCatch({
      wav_file <- readWave(file_path, from = read_from, to = read_to, units = "seconds")
    }, error = function(e) {
      stop(paste("Erreur lecture WAV:", e$message))
    })
    
  } else {
    # --- CAS AUTRES FORMATS (MP3, FLAC, etc.) ---
    temp_wav_path <- tempfile(fileext = ".wav")
    success <- tryCatch({
      av::av_audio_convert(file_path, temp_wav_path)
      TRUE
    }, error = function(e) {
      message("Erreur conversion: ", e$message)
      FALSE
    })
    
    if (!success || !file.exists(temp_wav_path)) {
      if (file.exists(temp_wav_path)) unlink(temp_wav_path)
      stop("Impossible de convertir le fichier audio.")
    }
    
    file_info <- get_audio_info(temp_wav_path)
    
    if (is.null(file_info)) {
      unlink(temp_wav_path)
      stop("Impossible de lire le fichier converti.")
    }
    
    Samp_rate <- file_info$samp.rate
    bit_depth_val <- file_info$bit
    full_duration <- file_info$duration
    
    read_to <- min(req_end, full_duration)
    read_from <- max(0, req_start)
    
    if (read_from >= read_to) {
      read_to <- min(full_duration, read_from + min_display_duration)
      read_from <- max(0, read_to - min_display_duration)
    }
    
    tryCatch({
      wav_file <- readWave(temp_wav_path, from = read_from, to = read_to, units = "seconds")
      unlink(temp_wav_path) # Nettoyage immédiat
    }, error = function(e) {
      unlink(temp_wav_path)
      stop(paste("Erreur lecture fichier converti:", e$message))
    })
  }
  
  # --- 4. VÉRIFICATIONS FINALES ET TRAITEMENT ---
  
  if (is.null(wav_file)) {
    stop("Échec du chargement du fichier audio.")
  }
  
  # Recalcul de la durée réelle du segment chargé (au cas où)
  # Pour un objet Wave, la durée est length(left) / samp.rate
  if (inherits(wav_file, "Wave")) {
    wav_duration <- length(wav_file@left) / wav_file@samp.rate
    Samp_rate <- wav_file@samp.rate # On récupère la valeur sûre de l'objet chargé
    bit_depth_val <- wav_file@bit
  } else {
    # Fallback très rare si l'objet n'est pas reconnu comme Wave
    stop("Le fichier chargé n'est pas un objet Wave valide.")
  }
  
  if (wav_duration < 0.5) {
    stop("Le segment audio est trop court pour être affiché (< 0.5s).")
  }
  
  # Découpe finale (ici, c'est souvent redondant car readWave a déjà coupé, mais sécurise)
  wav_segment <- cutw(wav_file, from = 0, to = wav_duration, f = Samp_rate, output = "Wave")
  
  if (length(wav_segment@left) == 0) {
    stop("Segment audio vide.")
  }
  
  # --- 5. TRAITEMENT DU SIGNAL (OISEAUX) ---
  
  # A. Réduction de bruit
  wav_clean <- tryCatch({
    rmnoise(wav_segment, f = Samp_rate, p = 0.85)
  }, error = function(e) {
    wav_segment@left
  })
  
  # B. Filtre Passe-Haut
  hp_freq <- max(100, highpass)
  wav_filtered <- tryCatch({
    fir(wav_clean, from = hp_freq, f = Samp_rate)
  }, error = function(e) {
    wav_clean
  })
  
  # Reconstruction objet Wave
  if (!is.numeric(wav_filtered)) wav_filtered <- as.numeric(wav_filtered)
  wav_final <- Wave(left = wav_filtered, samp.rate = Samp_rate, bit = bit_depth_val)


  
  # --- 6. NORMALISATION ---
  max_possible <- 2^(wav_final@bit - 1) - 1
  data_float <- wav_final@left / max_possible
  
  max_val <- max(abs(data_float), na.rm = TRUE)
  if (is.na(max_val) || max_val == 0) max_val <- 1
  
  target_level <- 10^(-6/20) 
  gain <- target_level / max_val
  wav_final@left <- data_float * gain
  
  if (length(wav_final@left) < 100 || all(is.na(wav_final@left))) {
    stop("Signal audio invalide après traitement.")
  }
  
  return(list(wav_final, Samp_rate, wav_file, wav_final@bit))

}

