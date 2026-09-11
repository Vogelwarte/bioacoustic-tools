
#' Get Root Directories Based on Operating System


get_roots <- function() {
  computer <- Sys.info()[["sysname"]]
  
  # Fonction helper pour extraire le nom propre d'un chemin
  # Ex: "/Volumes/MonDisque" -> "MonDisque"
  get_volume_name <- function(path) {
    basename(path)
  }
  
  if (computer == "Darwin") {
    # MacOS
    volumes <- list.dirs("/Volumes", full.names = TRUE, recursive = FALSE)
    
    # Exclusions
    exclusions <- c("/Volumes", "/Volumes/.timemachine", "/Volumes/com.apple.TimeMachine.localsnapshots")
    volumes <- setdiff(volumes, exclusions)
    
    # Création des noms propres
    volume_names <- basename(volumes) # Ex: "MonDisque"
    
    # On combine Home et les Volumes
    # Home sera nommé "Home Directory" ou "~"
    all_paths <- c("~", volumes)
    all_names <- c("Home Directory", volume_names)
    
    # Nettoyage des noms vides (cas rare)
    all_names[all_names == ""] <- "Unknown Volume"
    
    return(setNames(all_paths, all_names))
    
  } else if (computer == "Windows") {
    # Windows : Détection native
    tryCatch({
      # --- MÉTHODE ROBUSTE : TEST DIRECT DES LETTRES DE LECTEURS ---
      # On teste toutes les lettres de A à Z. C'est lent (quelques ms) mais 100% fiable.
      letters_avail <- LETTERS[1:26]
      found_paths <- character(0)
      found_names <- character(0)
      
      for (l in letters_avail) {
        drive_path <- paste0(l, ":/")
        
        # On vérifie si le dossier racine existe et est accessible
        if (dir.exists(drive_path)) {
          found_paths <- c(found_paths, drive_path)
          
          # Tentative de récupération du nom du volume (Label)
          # On essaie d'abord une méthode simple via list.files (parfois le nom apparait)
          # Sinon on met un nom par défaut "Disque X"
          vol_name <- paste0("Disque ", l)
          
          found_names <- c(found_names, vol_name)
        }
      }
      
      if (length(found_paths) > 0) {
        # Création de la liste nommée
        drives_list <- setNames(found_paths, found_names)
        
        # Nettoyage des noms vides ou doublons
        drives_list <- drives_list[nzchar(names(drives_list))]
        drives_list <- drives_list[!duplicated(drives_list)]
        
        return(drives_list)
        
      } else {
        stop("Aucun lecteur détecté de A: à Z:.")
      }
      
    }, error = function(e) {
      warning("Erreur détection disques: ", e$message, ". Retour secours sur C:.")
      return(c("Disque Local" = "C:/"))
    })
    
  } else if (computer == "Linux") {
    # Linux
    volumes_media <- list.dirs("/media", full.names = TRUE, recursive = FALSE)
    volumes_mnt   <- list.dirs("/mnt", full.names = TRUE, recursive = FALSE)
    volumes <- unique(c(volumes_media, volumes_mnt))
    
    # Exclusions
    volumes <- volumes[!volumes %in% c("/media", "/mnt")]
    
    # Noms propres
    volume_names <- basename(volumes)
    volume_names[volume_names == ""] <- "Unknown Volume"
    
    all_paths <- c("~", volumes)
    all_names <- c("Home Directory", volume_names)
    
    return(setNames(all_paths, all_names))
    
  } else {
    return(c(Home = "~"))
  }
}
