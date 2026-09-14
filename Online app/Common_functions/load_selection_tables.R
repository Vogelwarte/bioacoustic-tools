# 
# 
# 
# load_selection_tables <- function(dir1, dir2 = NULL, compiled = FALSE, utc_tz = "UTC") {
dir1 <- "/Volumes/Untitled/Results/"
#   # Liste des fichiers
#   pattern <- if (compiled) "*BirdNET_SelectionTable.txt" else "*.selection.table.txt"
#   filelist <- dir_ls(path = dir1, glob = pattern, recurse = TRUE)
#   
#   if (length(filelist) == 0) {
#     return(NULL)  # ou stop("Aucun fichier trouvé")
#   }
#   
#   # Parallélisation
#   nb_cores <- detectCores() - 1
#   cl <- makeCluster(nb_cores)
#   on.exit(stopCluster(cl), add = TRUE)
#   
#   DT_List <- parLapply(cl, filelist, function(f) {
#     DT_tmp <- data.table::fread(f, check.names = TRUE)
#     if (!"Begin.Path" %in% names(DT_tmp) & !"FileName" %in% names(DT_tmp)) {
#       DT_tmp[, file_path := f]
#     }
#     DT_tmp
#   })
#   
#   DT <- rbindlist(DT_List)
#   
#   # Encodage UTF-8
#   Encoding(DT$Common.Name) <- Encoding(enc2utf8(DT$Common.Name))
#   
#   # Gestion des chemins
#   if ("Begin.Path" %in% names(DT)) {
#     DT$Begin_Path <- DT$Begin.Path
#     DT$Begin_Path <- gsub("\\\\", "/", DT$Begin_Path)
#   } else {
#     if (!"FileName" %in% names(DT)) {
#       DT$FileName <- DT$file_path
#     }
#     DT$Begin_Path <- gsub(".BirdNET.selection.table.txt", ".WAV", DT$FileName, ignore.case = TRUE)
#   }
#   
#   # Extraction du nom de fichier sans extension
#   DT$file_name_basic <- str_remove(basename(DT$Begin_Path), "\\.(wav|WAV|mp3|Mp3|MP3|flac|FLAC)$")
#   
#   # Datetime
#   DT$DateTime_2 <- as.POSIXct(str_sub(DT$file_name_basic, start = -15), 
#                               format = "%Y%m%d_%H%M%S", tz = utc_tz)
#   
#   # Heure, date, minute
#   if (compiled) {
#     DT$Hour <- hour(DT$DateTime_2 + DT$File.Offset..s.)
#     DT$Date <- as.Date(DT$DateTime_2 + DT$File.Offset..s.)
#     DT$Min <- minute(DT$DateTime_2 + DT$File.Offset..s.)
#     DT$Start_segment <- DT$DateTime_2 + DT$File.Offset..s.
#     DT$Stop_segment <- DT$DateTime_2 + DT$File.Offset..s. + (DT$End.Time..s. - DT$Begin.Time..s.)
#   } else {
#     DT$Hour <- hour(DT$DateTime_2 + DT$Begin.Time..s.)
#     DT$Date <- as.Date(DT$DateTime_2 + DT$Begin.Time..s.)
#     DT$Min <- minute(DT$DateTime_2 + DT$Begin.Time..s.)
#     DT$Start_segment <- DT$DateTime_2 + DT$Begin.Time..s.
#     DT$Stop_segment <- DT$DateTime_2 + DT$End.Time..s.
#   }
#   
#   # Exclusion des "nocall"
#   DT <- DT[Common.Name != "nocall", ]
#   
#   # Association des fichiers audio (si dir2 fourni)
#   if (!is.null(dir2)) {
#     audio_files <- data.frame(
#       path = dir_ls(path = dir2, regexp = "\\.(wav|WAV|mp3|Mp3|MP3|flac|FLAC)$", 
#                     recurse = TRUE, ignore.case = TRUE)
#     )
#     audio_files$file_name_basic <- str_remove(basename(audio_files$path), 
#                                               "\\.(wav|WAV|mp3|Mp3|MP3|flac|FLAC)$")
#     
#     # Jointure
#     DT$True_Location <- audio_files$path[match(DT$file_name_basic, audio_files$file_name_basic)]
#   } else {
#     DT$True_Location <- NA_character_
#   }
#   
#   return(DT)
# }


# Pour la parallélisation moderne (souvent plus rapide et gère mieux la mémoire)
# Si non installé : install.packages("future.apply")


load_selection_tables <- function(dir1, dir2 = NULL, compiled = FALSE, utc_tz = "UTC") {
  
  # 1. Configuration initiale
  pattern <- if (compiled) "*BirdNET_SelectionTable.txt" else "*.selection.table.txt"
  filelist <- dir_ls(path = dir1, glob = pattern, recurse = TRUE)
  
  if (length(filelist) == 0) {
    message("Aucun fichier trouvé.")
    return(NULL)
  }
  
  message(sprintf("Traitement de %d fichiers...", length(filelist)))
  
  # 2. Stratégie de parallélisation
  # On planifie le backend. 'multisession' est souvent plus efficace que 'fork' sur Windows/Mac
  plan(multisession, workers = max(1, detectCores() - 1))
  
  # Fonction de lecture individuelle optimisée
  read_one_file <- function(f) {
    # fread est déjà très rapide. check.names=TRUE est parfois lent, essayez FALSE si vos noms sont propres
    # na.strings = "" évite les warnings inutiles
    DT_tmp <- fread(f, check.names = TRUE, na.strings = "", encoding = "UTF-8") 
    
    if (nrow(DT_tmp) == 0) return(NULL)
    
    # Ajout du chemin si colonnes manquantes (logique conservée)
    if (!"Begin.Path" %in% names(DT_tmp) && !"FileName" %in% names(DT_tmp)) {
      DT_tmp[, file_path := f]
    }
    return(DT_tmp)
  }
  
  # Exécution parallèle
  # future_lapply retourne une liste, gère automatiquement la progression si progressr est utilisé
  DT_List <- future_lapply(filelist, read_one_file)
  
  # Nettoyage du plan parallèle
  plan(sequential)
  
  # Filtrer les NULL (fichiers vides) et combiner
  DT_List <- DT_List[!sapply(DT_List, is.null)]
  
  if (length(DT_List) == 0) return(NULL)
  
  # rbindlist est extrêmement rapide
  DT <- rbindlist(DT_List, fill = TRUE, use.names = TRUE)
  
  # 3. Nettoyage et Transformation (Par référence avec := pour éviter les copies)
  
  # Gestion des chemins (Logique simplifiée)
  if ("Begin.Path" %in% names(DT)) {
    DT[, Begin_Path := gsub("\\\\", "/", Begin.Path)]
  } else {
    if (!"FileName" %in% names(DT)) {
      # Fallback si aucune colonne de nom n'existe
      if ("file_path" %in% names(DT)) {
        DT[, FileName := file_path]
      } else {
        stop("Impossible de déterminer le nom du fichier source.")
      }
    }
    # Correction du pattern de remplacement (plus robuste)
    DT[, Begin_Path := sub("\\.BirdNET\\.selection\\.table\\.txt$", ".WAV", FileName, ignore.case = TRUE)]
    # Si le nom ne correspondait pas, Begin_Path sera identique à FileName, ce qui peut être voulu ou non.
  }
  
  # Extraction nom de base (optimisé avec sub au lieu de str_remove pour la vitesse sur gros vecteur)
  DT[, file_name_basic := sub("\\.(wav|WAV|mp3|Mp3|MP3|flac|FLAC)$", "", basename(Begin_Path))]
  
  # 4. Gestion des Dates (Vectorisé et calcul unique)
  # Extraction de la date depuis le nom de fichier
  # On suppose que les 15 derniers caractères sont toujours YYYYMMDD_HHMMSS
  # Si la longueur varie, il faut une regex plus complexe, mais sub est plus rapide que str_sub ici
  date_strings <- substr(DT$file_name_basic, nchar(DT$file_name_basic) - 14, nchar(DT$file_name_basic))
  
  DT[, DateTime_2 := as.POSIXct(date_strings, format = "%Y%m%d_%H%M%S", tz = utc_tz)]
  
  # Calcul du décalage une seule fois
  if (compiled) {
    DT[, time_offset := File.Offset..s.]
  } else {
    # Vérifier si la colonne existe pour éviter l'erreur
    if ("Begin.Time..s." %in% names(DT)) {
      DT[, time_offset := Begin.Time..s.]
    } else {
      DT[, time_offset := 0] # Fallback si colonne manquante
    }
  }
  
  # Calcul absolu une seule fois
  DT[, Start_segment := DateTime_2 + time_offset]
  
  # Calcul des colonnes dérivées à partir de Start_segment (beaucoup plus rapide)
  DT[, `:=`(
    Date = as.Date(Start_segment),
    Hour = hour(Start_segment),
    Min = minute(Start_segment)
  )]
  
  # Calcul Stop_segment
  if (compiled) {
    DT[, Stop_segment := Start_segment + (End.Time..s. - Begin.Time..s.)]
  } else {
    if ("End.Time..s." %in% names(DT)) {
      DT[, Stop_segment := DateTime_2 + End.Time..s.]
    } else {
      DT[, Stop_segment := Start_segment] # Fallback
    }
  }
  
  # Suppression des colonnes temporaires pour libérer de la mémoire immédiatement
  #DT[, c("time_offset", "date_strings") := NULL]
  
  # Exclusion nocall (très rapide avec data.table)
  DT <- DT[Common.Name != "nocall"]
  
  # 5. Jointure Audio (Data Table Join binaire - Ultra rapide)
  if (!is.null(dir2)) {
    message("Association des fichiers audio...")
    audio_files <- data.table(
      path = dir_ls(path = dir2, regexp = "\\.(wav|WAV|mp3|Mp3|MP3|flac|FLAC)$", 
                    recurse = TRUE, ignore.case = TRUE)
    )
    audio_files[, file_name_basic := sub("\\.(wav|WAV|mp3|Mp3|MP3|flac|FLAC)$", "", basename(path))]
    
    # Définir les clés pour la jointure
    setkey(DT, file_name_basic)
    setkey(audio_files, file_name_basic)
    
    # Jointure directe (mise à jour par référence)
    DT[audio_files, True_Location := i.path]
  } else {
    DT[, True_Location := NA_character_]
  }
  

  
  return(DT)
}


