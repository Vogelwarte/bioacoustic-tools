
# Function to remove accents, punctuation and special characters to rename files for export
clean_text <- function(vecteur) {
  # Remove accents
  vecteur_sans_accents <- stringi::stri_trans_general(vecteur, "Latin-ASCII")
  
  # Remove punctuation and special characters
  vecteur_sans_speciaux <- stringr::str_replace_all(vecteur_sans_accents, "[^[:alnum:] ]", "")
  
  return(vecteur_sans_speciaux)
}

