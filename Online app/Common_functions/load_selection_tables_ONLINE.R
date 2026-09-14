load_selection_tables <- function(dir1, dir2 = NULL, compiled = TRUE,
                                  device_tz = "Europe/Zurich",  # Fuseau de l'horloge interne (ex: Alaska)
                                  deployment_tz = "Europe/Zurich")  # Fuseau du lieu réel (ex: France)
{
  # 1. Configuration initiale

  pattern <- if (compiled) "*BirdNET_SelectionTable.txt" else "*.selection.table.txt"
  filelist <- dir_ls(path = dir1, glob = pattern, recurse = TRUE)

  if (length(filelist) == 0) {
    message("Aucun fichier trouvé.")
    return(NULL)
  }

  message(sprintf("Traitement de %d fichiers...", length(filelist)))
  message(sprintf(" -> Heure interne (Device): %s", device_tz))
  message(sprintf(" -> Heure réelle (Deployment): %s", deployment_tz))

  # 2. Lecture séquentielle (sans parallélisation)
  read_one_file <- function(f) {
    DT_tmp <- fread(f, check.names = TRUE, na.strings = "", encoding = "UTF-8")
    if (nrow(DT_tmp) == 0) return(NULL)
    if (!"Begin.Path" %in% names(DT_tmp) && !"FileName" %in% names(DT_tmp)) {
      DT_tmp[, file_path := f]
    }
    return(DT_tmp)
  }

  DT_List <- lapply(filelist, read_one_file)
  DT_List <- DT_List[!sapply(DT_List, is.null)]
  if (length(DT_List) == 0) return(NULL)
  DT <- rbindlist(DT_List, fill = TRUE, use.names = TRUE)

  # 3. Nettoyage et Transformation
  if ("Begin.Path" %in% names(DT)) {
    DT[, Begin_Path := gsub("\\\\", "/", Begin.Path)]
  } else {
    if (!"FileName" %in% names(DT)) {
      if ("file_path" %in% names(DT)) DT[, FileName := file_path]
      else stop("Impossible de déterminer le nom du fichier source.")
    }
    DT[, Begin_Path := sub("\\.BirdNET\\.selection\\.table\\.txt$", ".WAV", FileName, ignore.case = TRUE)]
  }

  DT[, file_name_basic := sub("\\.(wav|WAV|mp3|Mp3|MP3|flac|FLAC)$", "", basename(Begin_Path))]

  # 4. Gestion des Dates - ÉTAPE 1 : LECTURE DE L'HEURE "DISPLAY" (DEVICE)
  # Extraction de la chaîne date_heure (les 15 derniers caractères)
  DT[, date_strings := substr(file_name_basic, pmax(1, nchar(file_name_basic) - 14), nchar(file_name_basic))]

  # Création de DateTime_Display dans le fuseau de l'appareil (ex: Alaska)
  DT[, DateTime_Display := as.POSIXct(date_strings, format = "%Y%m%d_%H%M%S", tz = device_tz)]
  #attr(DT$DateTime_Display, "tzone") <- device_tz # Verrouillage

  # Colonnes dérivées "Display" (pour référence/logs)
  DT[, `:=`(
    Date_Display = as.Date(DateTime_Display),
    Hour_Display = hour(DateTime_Display),
    Min_Display = minute(DateTime_Display),
    Time_Display = sprintf("%02d:%02d", hour(DateTime_Display), minute(DateTime_Display))
  )]

  # 5. Gestion des Dates - ÉTAPE 2 : CONVERSION VERS L'HEURE "RÉELLE" (DEPLOYMENT)
  if (!is.null(deployment_tz) && deployment_tz != device_tz) {
    # Conversion via UTC pour garantir la justesse de l'instant
    # 1. On passe en UTC (instant universel)
    utc_strings <- format(DT$DateTime_Display, tz = "UTC", usetz = FALSE)

    # 2. On recrée l'objet dans le fuseau de déploiement (ex: France)
    DT[, DateTime_Real := as.POSIXct(utc_strings, format = "%Y-%m-%d %H:%M:%S", tz = deployment_tz)]
    attr(DT$DateTime_Real, "tzone") <- deployment_tz

    message("Conversion des fuseaux horaires effectuée.")
  } else {
    # Si mêmes fuseaux, on duplique
    DT[, DateTime_Real := DateTime_Display]
    message("Pas de conversion de fuseau (identiques).")
  }


  # 6. Calcul du segment temporel absolu (Start_segment)
  # On utilise l'heure RÉELLE pour les calculs scientifiques (soleil, etc.)
  if (compiled) {
    DT[, time_offset := File.Offset..s.]
  } else {
    if ("Begin.Time..s." %in% names(DT)) {
      DT[, time_offset := Begin.Time..s.]
    } else {
      DT[, time_offset := 0]
    }
  }

  # Calcul de Start_segment basé sur DateTime_Real
  DT[, Start_segment := DateTime_Real + time_offset]

  # Calcul Stop_segment
  if (compiled) {
    if ("End.Time..s." %in% names(DT)) {
      DT[, Stop_segment := Start_segment + (End.Time..s. - File.Offset..s.)]
    } else {
      DT[, Stop_segment := Start_segment]
    }
  } else {
    if ("End.Time..s." %in% names(DT)) {
      DT[, Stop_segment := DateTime_Real + End.Time..s.]
    } else {
      DT[, Stop_segment := Start_segment]
    }
  }

  DT[, `:=`(
    Date = as.Date(Start_segment),
    Hour = hour(Start_segment),
    Min = minute(Start_segment),
    Time = sprintf("%02d:%02d", hour(Start_segment), minute(Start_segment)),
    Hour_Decimal_Real = hour(Start_segment) + (minute(Start_segment) / 60)
  )]


  # Nettoyage colonne temporaire
  DT[, time_offset := NULL]

  # Exclusion nocall
  if ("Common.Name" %in% names(DT)) {
    DT <- DT[Common.Name != "nocall"]
  }

  # 7. Jointure Audio
  if (!is.null(dir2)) {
    message("Association des fichiers audio...")
    audio_files <- data.table(
      path = dir_ls(path = dir2, regexp = "\\.(wav|WAV|mp3|Mp3|MP3|flac|FLAC)$", recurse = TRUE, ignore.case = TRUE)
    )
    audio_files[, file_name_basic := sub("\\.(wav|WAV|mp3|Mp3|MP3|flac|FLAC)$", "", basename(path))]

    setkey(DT, file_name_basic)
    setkey(audio_files, file_name_basic)
    DT[audio_files, True_Location := i.path]
  } else {
    DT[, True_Location := NA_character_]
  }

  # 8. Parsing Robuste du Nom de l'Enregistreur
  # Extrait tout ce qui précède la date (YYYYMMDD) et l'heure (HHMMSS)
  DT[, recorder := sub("^(.*)_\\d{8}_\\d{6}($|\\..*)$", "\\1", file_name_basic)]

  # Fallback si le regex échoue (nom de fichier atypique) : prend le premier élément
  DT[recorder == file_name_basic, recorder := sapply(strsplit(file_name_basic, "_"), `[`, 1)]

  return(DT)
}
