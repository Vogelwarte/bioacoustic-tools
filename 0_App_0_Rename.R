# ---- Install if needed ----
# install.packages(c("shiny", "shinyFiles", "leaflet"))

library(shiny)
library(shinyFiles)
library(leaflet)
library(bslib)

ui <- fluidPage(
  # --- En-tête personnalisé avec Logo ---
  div(
    style = "display: flex; justify-content: space-between; align-items: center; padding: 10px 20px; background-color: #0f0f0f; border-bottom: 1px solid #333; margin-bottom: 20px;",
    div(
      style = "font-size: 24px; font-weight: bold; color: #fff;",
      "BirdNET-ResChecker - Rename Files with Coordinates and Recorder Prefix"
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
  
  #titlePanel("Rename Files with Coordinates and Recorder Prefix"),
  
  sidebarLayout(
    sidebarPanel(
      shinyDirButton("folder", "Choose Folder", "Select a folder"),
      textInput("prefix", "Prefix to add:", value = ""),
      actionButton("preview", "Preview Changes"),
      textInput("prefix", "Recommended synthax: Lat_Long_Recorder_Date_Time.wav"),
      actionButton("rename", "Apply Renaming"),
      tags$hr(),
      h4("Selected Coordinates"),
      verbatimTextOutput("coords_text"),
      tags$hr(),
      verbatimTextOutput("selected_folder"),
      verbatimTextOutput("status")
    ),
    
    mainPanel(
      tabsetPanel(type = "tabs",
                  
                  tabPanel("README",
                           h3("Documentation"),
                           p("This script enables batch renaming of files by adding a location-based prefix. It processes and renames all files within a selected directory."),
                           p("The prefix can be entered manually or generated automatically by clicking on the map 🌍."),
                           p("⚠️ Important ⚠️ : Files names should follow this format : LatXX.XXXX_LongYY.YYYY_Recorder_Date_Time.wav"),
                           
                           
                           h4("Usage Instructions"),
                           tags$ol(
                             tags$li("Click 'Choose Folder' to select the directory containing the files."),
                             tags$li("Generate a prefix by either:"),
                             tags$ul("Clicking on the map to automatically retrieve latitude/longitude, or"),
                             tags$ul("Entering the prefix manually."),
                             tags$li("⚠️ If entering manually, ensure the correct syntax: LatXX.XXXX_LongYY.YYYY_Recorder_Date_Time.wav"),
                             tags$li("Click 'Preview Changes' to visualize the updated file names."),
                             tags$li("Click 'Apply Renaming' to perform the renaming operation. \n ⚠️ Warning: This action is irreversible. Original file names cannot be restored once the process is executed.")
                           ),
                           h4("Dependencies"),
                           p("Required R packages: shiny, shinyFiles, leaflet, bslib"),
                           #p("To install, run:"),
                           #verbatimTextOutput("install_cmd")
                  ),
                  tabPanel("Map & Rename",
                           h4("Select location on the map to auto-generate prefix"),
                           leafletOutput("map", height = "400px"),
                           tags$hr(),
                           h4("Preview of Changes"),
                           tableOutput("preview_table")
                  ),
                  nav_panel("Reference",
                            layout_columns(
                              col_widths = c(2, 8, 2),
                              textOutput("contributions"),
                              textOutput("license")
                            )
      ),
      )
    )
  )
)

server <- function(input, output, session) {
  
  # --- Textes statiques ---
  output$contributions <- renderText("Contributions: Amandine Serrurier, Jean-Nicolas Pradervand & Christophe Sahli\nSwiss Ornithological Institute")
  output$license <- renderText("MIT License © 2026 Jean-Nicolas Pradervand")
  
  # --- Detect platform and set roots dynamically ---
  if (.Platform$OS.type == "windows") {
    
    # All available drives
    volumes <- getVolumes()()
    roots <- c(Home = normalizePath("~"), volumes)
    
  } else if (Sys.info()[["sysname"]] == "Darwin") {
    
    # macOS: external drives live in /Volumes
    roots <- c(
      Home = normalizePath("~"),
      Volumes = "/Volumes",
      Root = "/"
    )
    
  } else {
    
    # Linux
    roots <- c(
      Home = normalizePath("~"),
      Media = "/media",
      Mnt = "/mnt",
      Root = "/"
    )
  }
  
  # --- Folder selection setup ---
  shinyDirChoose(input, "folder", roots = roots, session = session)
  
  folder_path <- reactive({
    req(input$folder)
    parseDirPath(roots, input$folder)
  })
  
  output$selected_folder <- renderText({
    req(folder_path())
    paste("Selected folder:", folder_path())
  })
  
  # --- Map setup ---
  output$map <- renderLeaflet({
    leaflet() |>
      addTiles() |>
      setView(lng = 10, lat = 50, zoom = 3)
  })
  
  # --- When user clicks on map ---
  observeEvent(input$map_click, {
    click <- input$map_click
    lat <- round(click$lat, 5)
    lon <- round(click$lng, 5)
    prefix <- paste0("Lat", lat, "_Long", lon, "_")
    
    updateTextInput(session, "prefix", value = prefix)
    output$coords_text <- renderText({
      paste("Latitude:", lat, "\nLongitude:", lon)
    })
  })
  
  # --- Preview changes ---
  preview_data <- eventReactive(input$preview, {
    req(folder_path())
    req(input$prefix)
    files <- list.files(folder_path(), full.names = FALSE)
    if (length(files) == 0)
      return(data.frame(Message = "No files found in folder."))
    data.frame(
      Original = files,
      New = paste0(input$prefix, files),
      stringsAsFactors = FALSE
    )
  })
  
  output$preview_table <- renderTable({
    preview_data()
  })
  
  # --- Apply renaming ---
  observeEvent(input$rename, {
    req(folder_path())
    req(input$prefix)
    files <- list.files(folder_path(), full.names = TRUE)
    new_names <- file.path(folder_path(), paste0(input$prefix, basename(files)))
    success <- file.rename(files, new_names)
    
    if (all(success)) {
      output$status <- renderText("✅ All files renamed successfully.")
    } else {
      failed <- basename(files[!success])
      output$status <- renderText(
        paste("⚠️ Some files could not be renamed:", paste(failed, collapse = ", "))
      )
    }
  })
}

shinyApp(ui, server)
