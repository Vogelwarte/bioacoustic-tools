
#' Get Root Directories Based on Operating System


get_roots <- function() {
  computer <- Sys.info()[["sysname"]]
  
  if (computer == "Darwin") {
    # MacOS
    volumes <- list.dirs("/Volumes", full.names = TRUE, recursive = FALSE)
    # exclusions spécifiques à macOS
    exclusions <- c("/Volumes", 
                    "/Volumes/.timemachine", 
                    "/Volumes/com.apple.TimeMachine.localsnapshots")
    volumes <- setdiff(volumes, exclusions)
    volumes <- c("~", volumes)
    return(setNames(volumes, volumes))
    
  } else if (computer == "Windows") {
    # Windows (liste personnalisable)
    drives <- c("C:/", "D:/", "E:/", "F:/","G:/","H:/", "I:/", "J:/" , "M:/", "V:/")
    return(setNames(drives, drives))
    
  } else if (computer == "Linux") {
    # Linux
    volumes_media <- list.dirs("/media", full.names = TRUE, recursive = FALSE)
    volumes_mnt   <- list.dirs("/mnt", full.names = TRUE, recursive = FALSE)
    volumes <- unique(c(volumes_media, volumes_mnt))
    volumes <- volumes[!volumes %in% c("/media", "/mnt")]
    return(setNames(volumes, volumes))
    
  } else {
    warning("Unsupported operating system.")
    return(NA)
  }
}


