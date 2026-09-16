#' Charger et nettoyer les données brutes BirdNET
#'
#' @param files_obj L'objet input$... provenant de fileInput (ex: input$txtfiles)
#' @return Un data.frame tidy (tibble) ou NULL si aucun fichier
#' @importFrom purrr map2_dfr map_chr map
#' @importFrom readr read_delim cols parse_number
#' @importFrom dplyr mutate bind_rows
#' @importFrom stringr str_trim str_to_lower str_extract str_split str_remove str_detect
#' @importFrom janitor clean_names
load_raw_birdnet_data <- function(files_obj) {
  #files_obj <- DT
  # Vérification de sécurité si aucun fichier
  if (is.null(files_obj) || (is.data.frame(files_obj) && nrow(files_obj) == 0)) {
    return(NULL)
  }
  
  
  data <- files_obj %>%
    #janitor::clean_names() %>%
    dplyr::mutate(
      common_name_original = Common.Name,
      common_name = stringr::str_trim(stringr::str_to_lower(Common.Name)),
      
      # Détermination du nom de fichier à utiliser
      # Si la colonne 'begin_path' existe, on prend son basename, sinon on prend 'source_file'
      filename = if ("file_name_basic" %in% names(.)) {
        file_name_basic
      } else {
        basename(Begin.Path)
      },
      
      filename_clean = stringr::str_remove(
        filename,
        "^Lat-?\\d+\\.\\d+_Long-?\\d+\\.\\d+_"
      ),
      
      # Extraction Lat/Long
      lat = readr::parse_number(stringr::str_extract(filename, "Lat-?\\d+\\.\\d+")),
      long = readr::parse_number(stringr::str_extract(filename, "Long-?\\d+\\.\\d+")),
      
      # Découpage du nom de fichier
      parts = stringr::str_split(filename_clean, "_"),
      
      # Nettoyage de la liste des parties
      parts_clean = purrr::map(parts, function(x) {
        if (length(x) >= 3 && 
            stringr::str_detect(x[1], "^Lat") && 
            stringr::str_detect(x[2], "^Long")) {
          return(x[-c(1, 2)])
        } else {
          return(x)
        }
      }),
      
      # Extraction des métadonnées
      recorder = as.factor(purrr::map_chr(parts_clean, function(x) {
        if (length(x) >= 1) return(x[1]) else return(NA_character_)
      })),
      
      date_str = purrr::map_chr(parts_clean, function(x) {
        if (length(x) >= 2) return(x[2]) else return(NA_character_)
      }),
      
      time_str = purrr::map_chr(parts_clean, function(x) {
        if (length(x) >= 3) return(stringr::str_remove(x[3], "\\.wav$|\\.txt$")) else return(NA_character_)
      }),
      
      # Conversion Date/Heure
      datetime = as.POSIXct(
        paste0(date_str, time_str),
        format = "%Y%m%d%H%M%S",
        tz = "UTC"
      ),
      
      date = as.Date(date_str, "%Y%m%d")
    )
  print(unique(data$recorder)) # test for recorder names
  return(data)
}