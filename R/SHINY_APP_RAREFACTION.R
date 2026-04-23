#=============================================================================#
# SHINY_APP_RAREFACTION.R
# Acoustic Sampling Rarefaction — Shiny Application
#
# Tabs:
#   1. Data import      – load & filter BirdNET output files
#   2. Duty cycle       – rarefaction across duty-cycle periods
#   3. Time window      – rarefaction across recording window start / duration
#   4. Multi optimum    – combined duty-cycle × time-window Pareto optimum
#   5. Target species   – temporal activity heatmap for focal species
#   User Guide          – in-app documentation
# =============================================================================#

# Libraries ----

library(shiny)
library(shinythemes)
library(data.table)
library(lubridate)
library(ggplot2)
library(suncalc)
library(dplyr)
library(DT)
library(scico)
library(plotly)
library(stringr)


# Shiny options ----

options(
  shiny.maxRequestSize  = 8 * 1024^3,   # 8 GB upload limit
  future.globals.maxSize = 8 * 1024^3   # 8 GB limit for parallel processing
)

# Load functions ----
source("functions.R")


# Recorder battery specifications----

battery_specs <- data.frame(
  recorder     = c("SMMicro 2", "SMMini 2", "SMMini 2", "SM4", "SM5", "AudioMoth -v1.2.0"),
  nb_battery   = c(4, 8, 6, 4, 11, 3),
  battery_type = c("AA", "AA", "lithium", "4D alkaline", "lithium", "AA"),
  duration_h   = c(280, 530, 1330, 310, 1710, 189)
)
battery_specs$code <- paste0(
  gsub(" ", "", battery_specs$recorder), "_",
  battery_specs$nb_battery,
  substr(battery_specs$battery_type, 1, 2)
)


# UI ----

ui <- fluidPage(
  theme = shinytheme("darkly"),
  
  ## Header ----
  div(
    style = "margin-top: 0px; padding-top: 0px;",
    fluidRow(
      style = "margin-top: 0px; padding-top: 0px;",
      column(
        width = 8,
        div(
          style = "margin-top: 5px; padding-top: 0px;",
          h2("Acoustic sampling rarefaction", style = "margin-top: 0px; margin-bottom: 0px;")
        )
      ),
      column(
        width = 4,
        div(
          style = "text-align: right; margin-top: 5px; padding-top: 10px;",
          div(
            style = "
              display: inline-block;
              background-color: white;
              padding: 6px 10px;
              border-radius: 8px;
              border: 1px solid #dddddd;
            ",
            tags$img(src = "logo.png", height = "70px")
          )
        )
      )
    )
  ),
  
  tabsetPanel(
    
    ## TAB 1 UI Data import----
    tabPanel(
      "1. Data import",
      sidebarLayout(
        sidebarPanel(
          
          h4("Import BirdNET data"),
          fileInput(
            inputId = "bn_file",
            label   = "Upload BirdNET file(s) (.csv, .txt, .zip)",
            accept  = c(".csv", ".txt", ".zip"),
            multiple = TRUE
          ),
          actionButton("load_data", "IMPORT AND PREVIEW"),
          
          hr(),
          
          h4("Path parsing (optional)"),
          div(
            style = "overflow-x: auto; max-width: 100%;",
            tableOutput("beginpath_preview")
          ),
          uiOutput("path_parser_ui"),
          br(),
          actionButton("apply_parsing", "Parse path"),
          
          hr(),
          h4("Date filter"),
          dateRangeInput(
            "date_subset",
            "Select time window"
          ),
          fluidRow(
            column(6, actionButton("apply_date_filter", "Filter dates")),
            column(6, actionButton("reset_date_filter", "Reset dates"))
          ),
          
          hr(),
          
          h4("Filter data"),
          numericInput("min_conf", "Min confidence", value = 0, min = 0, max = 1, step = 0.01),
          # dateRangeInput(
          #   "date_subset",
          #   "Select time window"
          #   # ,
          #   # start = "2024-06-01",
          #   # end   = "2024-06-09"
          # ),
          uiOutput("site_ui"),
          uiOutput("rec_ui"),
          uiOutput("year_ui"),
          uiOutput("species_ui"),
          actionButton("apply_filters", "Apply filters"),
          
          hr(),
          
          h4("Recorder coordinates"),
          radioButtons(
            "coord_mode",
            "Coordinate mode",
            choices = c(
              "Single coordinates (WGS84)"           = "single",
              "Add multiple coordinates (file upload)" = "multiple"
            ),
            selected = "single"
          ),
          selectizeInput(
            "timezone",
            "Timezone",
            choices  = OlsonNames(),
            selected = "CET",
            options  = list(placeholder = "Type to search timezone...", maxOptions = 1000)
          ),
          
          conditionalPanel(
            condition = "input.coord_mode == 'single'",
            numericInput("single_lat", "Latitude (WGS84)",  value = 46.19, step = 0.0001),
            numericInput("single_lon", "Longitude (WGS84)", value = 8.13,  step = 0.0001)
          ),
          
          conditionalPanel(
            condition = "input.coord_mode == 'multiple'",
            fileInput(
              "coord_file",
              "Upload coordinate table (csv, txt, xlsx)",
              accept = c(".csv", ".txt", ".xlsx")
            ),
            helpText("File must contain columns: recorder OR site + latitude + longitude (WGS84)")
          )
        ),
        mainPanel(
          h4("Preview"),
          tableOutput("table_preview")
        )
      )
    ),
    
    ## TAB2: UI Duty cycle ----
    tabPanel(
      "2. Duty cycle",
      fluidRow(
        column(
          width = 3,
          sidebarPanel(
            checkboxInput("record_24h", "Continuous 24h recording", FALSE),
            selectInput("start_event", "Start event", c("dawn", "sunrise", "sunset", "dusk")),
            selectInput("end_event",   "End event",   c("dawn", "sunrise", "sunset", "dusk")),
            numericInput("length_morning",  "Fixed duration (h)",    0, min = 0),
            numericInput("extend_before",   "Extend before (h)",     1),
            numericInput("extend_after",    "Extend after (h)",      0),
            
            checkboxInput("use_window2", "Add second daily window", FALSE),
            conditionalPanel(
              condition = "input.use_window2 == true",
              selectInput("start_event2",    "Start event (Window 2)", c("dawn", "sunrise", "sunset", "dusk")),
              selectInput("end_event2",      "End event (Window 2)",   c("dawn", "sunrise", "sunset", "dusk")),
              numericInput("length_morning2", "Fixed duration (h) — Window 2", 0, min = 0),
              numericInput("extend_before2",  "Extend before (h) — Window 2", 0),
              numericInput("extend_after2",   "Extend after (h) — Window 2",  0)
            ),
            
            numericInput("period_min",      "Min period (min)",   1),
            numericInput("period_max",      "Max period (min)",  30),
            numericInput("period_step",     "Step (min)",         5),
            numericInput("duty_duration",   "Duty duration (min)", 1),
            numericInput("B_temporal",      "Bootstrap",  value = 5, min = 1, step = 1),
            checkboxInput("show_summary",   "Show average + CI", TRUE),
            
            actionButton("run_temporal",          "Run analysis"),
            actionButton("return_raref",          "Return rarefaction plot",
                         style = "background-color: #728FCE; color: white; font-size: 12px; padding: 4px 8px; border: none;"),
            actionButton("return_heatmap_miss",   "Return missed species plot",
                         style = "background-color: #728FCE; color: white; font-size: 12px; padding: 4px 8px; border: none;"),
            actionButton("return_duty_data",      "Return duty cycle data",
                         style = "background-color: #728FCE; color: white; font-size: 12px; padding: 4px 8px; border: none;"),
            width = 12
          )
        ),
        column(
          width = 9,
          fluidRow(
            column(12, h4("Schedule preview")),
            column(6, plotOutput("sun_plot",  height = "150px")),
            column(6, plotOutput("duty_plot", height = "150px"))
          ),
          br(),
          plotlyOutput("plot_temporal",           height = "600px", width = "100%"),
          plotlyOutput("heatmap_missed_species",  height = "900px", width = "100%")
        )
      )
    ),
    
    ## TAB3: UI Time window ----
    tabPanel(
      "3. Time window",
      fluidRow(
        column(
          width = 3,
          sidebarPanel(
            selectInput(
              "grid_event",
              "Reference event",
              choices  = list(
                "Solar events" = c("dawn", "sunrise", "sunset", "dusk"),
                "Fixed hours"  = sprintf("%02d:00", 1:23)
              ),
              selected = "sunrise"
            ),
            numericInput("grid_start_min",  "Min offset from event (hours)",  value = -2,  step = 0.25),
            numericInput("grid_start_max",  "Max offset from event (hours)",  value =  4,  step = 0.25),
            numericInput("grid_start_step", "Start time step (hours)",        value =  0.5, step = 0.25),
            numericInput("grid_dur_min",    "Min duration (hours)",           value =  0.5, step = 0.25),
            numericInput("grid_dur_max",    "Max duration (hours)",           value =  6,  step = 0.25),
            numericInput("grid_dur_step",   "Duration step (hours)",          value =  0.5, step = 0.25),
            
            checkboxInput("use_window2_grid", "Add second daily window", FALSE),
            conditionalPanel(
              condition = "input.use_window2_grid == true",
              selectInput(
                "grid_event2",
                "Reference event (Window 2)",
                choices  = list(
                  "Solar events" = c("dawn", "sunrise", "sunset", "dusk"),
                  "Fixed hours"  = sprintf("%02d:00", 1:23)
                ),
                selected = "sunset"
              ),
              numericInput("grid_start_min2",  "Min offset (Window 2)",       value = -2,  step = 0.25),
              numericInput("grid_start_max2",  "Max offset (Window 2)",       value =  4,  step = 0.25),
              numericInput("grid_start_step2", "Start time step (Window 2)",  value =  0.5, step = 0.25),
              numericInput("grid_dur_min2",    "Min duration (Window 2)",     value =  0.5, step = 0.25),
              numericInput("grid_dur_max2",    "Max duration (Window 2)",     value =  6,  step = 0.25),
              numericInput("grid_dur_step2",   "Duration step (Window 2)",    value =  0.5, step = 0.25)
            ),
            
            actionButton("run_grid",          "Run analysis"),
            actionButton("return_window",     "Return day window plot",
                         style = "background-color: #728FCE; color: white; font-size: 12px; padding: 4px 8px; border: none;"),
            actionButton("return_window_data","Return day window data",
                         style = "background-color: #728FCE; color: white; font-size: 12px; padding: 4px 8px; border: none;"),
            width = 12,
            numericInput(
              "top_n_combined",
              "Keep top combined windows",
              value = 3,
              min   = 1,
              max   = 20,
              step  = 1
            ),
            checkboxInput("show_frontier", "Show optimal frontier", value = TRUE)
          )
        ),
        column(
          width = 9,
          h4("Design preview"),
          uiOutput("grid_summary"),
          DT::DTOutput("grid_preview"),
          br(),
          h4("Richness heatmap"),
          uiOutput("grid_plots_ui")
        )
      )
    ),
    
    ## TAB4: UI Multi optimum ----
    tabPanel(
      "4. Multi optimum",
      fluidRow(
        column(
          width = 12,
          helpText(
            "This analysis combines the optimized duty cycle and daily window.",
            "Please run the Temporal and Window analyses first."
          ),
          actionButton("run_multi",         "Run analysis"),
          actionButton("return_multi",      "Return multi-optimum plot",
                       style = "background-color: #728FCE; color: white; font-size: 12px; padding: 4px 8px; border: none;"),
          actionButton("return_multi_data", "Return multi optimum data",
                       style = "background-color: #728FCE; color: white; font-size: 12px; padding: 4px 8px; border: none;")
        )
      ),
      fluidRow(
        column(
          width = 12,
          plotlyOutput("plot_multi_optimum", height = "650px", width = "100%")
        )
      )
    ),
    
    ## Tab 5: UI Target species activity ----
    tabPanel(
      "5. Target species activity",
      sidebarLayout(
        sidebarPanel(
          width = 2,
          selectizeInput("target_species", "Target species", choices = NULL, multiple = TRUE),
          numericInput("target_period_min",      "Period min",             1,  min = 1),
          numericInput("target_period_max",      "Period max",            60,  min = 1),
          numericInput("target_period_step",     "Period step",            5,  min = 1),
          numericInput("target_duty_duration",   "Duty duration (min)",    1,  min = 1),
          numericInput("target_nb_duty",         "Number of duty per period", 1, min = 1),
          selectInput(
            "time_resolution",
            "Time resolution plot",
            choices  = c(5, 10, 15, 30, 45, 60),
            selected = 30
          ),
          actionButton("run_target",          "Run analysis"),
          actionButton("return_target",       "Return target species plot",
                       style = "background-color: #728FCE; color: white; font-size: 12px; padding: 4px 8px; border: none;"),
          actionButton("return_target_data",  "Return target species data",
                       style = "background-color: #728FCE; color: white; font-size: 12px; padding: 4px 8px; border: none;")
        ),
        mainPanel(
          width = 10,
          plotlyOutput("target_heatmap", height = "900px", width = "100%")
        )
      )
    ),
    
    # TAB User Guide ----
    tabPanel(
      "User Guide",
      fluidPage(
        tags$head(
          ##User guide CSS styles ----
          tags$style(HTML("
            .guide-main {
              background-color: #2b2b2b;
              border: 1px solid #444;
              border-radius: 12px;
              padding: 18px 22px;
              margin-bottom: 18px;
            }
            .guide-important {
              background-color: #3a2f1f;
              border-left: 4px solid #d9a441;
              border-radius: 8px;
              padding: 10px 14px;
              margin-top: 12px;
              margin-bottom: 12px;
            }
            .guide-highlight {
              background-color: #1f2d3a;
              border-left: 4px solid #4aa3df;
              border-radius: 8px;
              padding: 10px 14px;
              margin-top: 12px;
              margin-bottom: 12px;
            }
            details {
              background-color: #2b2b2b;
              border: 1px solid #444;
              border-radius: 10px;
              padding: 10px 14px;
              margin-bottom: 14px;
            }
            summary {
              font-weight: 600;
              font-size: 17px;
              cursor: pointer;
              outline: none;
            }
            .guide-subtitle {
              font-weight: 600;
              margin-top: 12px;
              margin-bottom: 4px;
            }
            .guide-img-box {
              text-align: center;
              padding: 6px;
            }
          "))
        ),
        ##User guide Application overview ----
        h2("Application overview"),
        
        div(
          class = "guide-main",
          fluidRow(
            column(
              width = 7,
              p(strong("Purpose")),
              tags$ul(
                tags$li("This application evaluates and optimises passive acoustic sampling designs from BirdNET detections."),
                tags$li("It uses an existing BirdNET dataset to simulate how different recording schedules would perform on the same data."),
                tags$li("It supports practical decisions on recording design while retaining a high proportion of species richness.")
              ),
              p(strong("What the app can test")),
              tags$ul(
                tags$li("Different duty-cycle settings."),
                tags$li("Different recording windows defined by start time and duration."),
                tags$li("Combinations of recording windows and duty-cycle settings."),
                tags$li("Temporal detection patterns of selected target species.")
              ),
              p(strong("Expected input")),
              tags$ul(
                tags$li("One BirdNET detection table, several BirdNET output files, or pre-merged files."),
                tags$li("The dataset should be representative of the study site and period of interest."),
                tags$li("In practice, a small pilot dataset with only a few recording days is sufficient.")
              ),
              div(
                class = "guide-highlight",
                p(strong("Suggested order of use")),
                tags$ol(
                  tags$li("Import and filter the BirdNET dataset. Define coordinates and timezone."),
                  tags$li("We suggest to start with a a small dataset of a subset of data as computations can take a long time, espcially if a large set of parameters are used."),
                  tags$li("Run the Duty cycle tab."),
                  tags$li("Run the Time window tab."),
                  tags$li("Run Multi optimum to compare complete sampling designs and identify the optimal frontier."),
                  tags$li("Use Target species activity to visualise temporal detection patterns of focal species.")
                )
              ),
              div(
                class = "guide-highlight",
                p(strong("Tips")),
                tags$ol(
                  tags$li("In each Tab, you can download data table associated with analysis and ggplot object. They will appear
                          in the global environement after closing the app."),
                  tags$li("In tabs refereing to recording effort, the effort of candidate schedule effort will be associated with battery life of some common recorders
                          for wildlife ans associated battery set. Informations about battery life are found on the website of each recorder.")
                )
              )
            ),
            column(
              width = 5,
              div(
                class = "guide-img-box",
                tags$img(src = "scheme_dark.png", style = "width: 80%; height: auto; border-radius: 10px;")
              )
            )
          )
        ),
        
        ##User guide by tab ----
        h2("Guide by tab", tags$span("(click to expand)", style = "font-size: 16px; font-weight: normal;")),
        ###User guide TAB1 ----
        tags$details(
          tags$summary("1. Data import"),
          p(class = "guide-subtitle", "Purpose of the tab"),
          tags$ul(
            tags$li("Import and prepare the BirdNET detection dataset used in all downstream analyses."),
            tags$li("BirdNET output is formatted to calculate the timestamp of each detection and to standardise column names."),
            tags$li("Assign site or recorder labels via path parsing."),
            tags$li("Filter by date, confidence score, species, recorder, or site."),
            tags$li("Assign coordinates and timezone to the recording location.")
          ),
          p(class = "guide-subtitle", "How to proceed"),
          tags$ol(
            tags$li("Browse one or several BirdNET files."),
            tags$li("Click ", strong("IMPORT AND PREVIEW"), "."),
            tags$li("Optional: use the path preview to identify where site and recorder names are stored in the file path. Assign the position of your information (e.g. V1 = site)."),
            tags$li("Select the corresponding path columns and click ", strong("Parse path"), "."),
            tags$li("Apply the desired filters."),
            tags$li("Assign one pair of coordinates for the whole dataset, or upload recorder-specific coordinates. Recorder-specific coordinates must be in a csv, txt, or xlsx file, and the location name must exactly match the name defined during path parsing."),
            tags$li("Select the timezone. Coordinates and timezone are used to calculate solar events (dawn, sunrise, sunset, dusk).")
          ),
          p(class = "guide-subtitle", "Accepted input"),
          tags$ul(
            tags$li("A single BirdNET detection table."),
            tags$li("Several BirdNET output files."),
            tags$li("Files already merged beforehand."),
            tags$li("Compressed .zip archives containing BirdNET selection tables.")
          ),
          div(
            class = "guide-important",
            p(strong("Important")),
            p("If path parsing is not applied, the same coordinates are used for the whole dataset. This is suitable when all recordings come from the same area, but it may reduce the accuracy of solar-event calculations if recorders are far apart."),
            p("If path parsing is applied, recorder- or site-specific coordinates must be provided in a multiple-coordinate file. Conversely, if path parsing is not applied, a single pair of coordinates must be used for the whole dataset.")
            )
        ),
        
        ###User guide TAB2 ----
        tags$details(
          tags$summary("2. Duty cycle"),
          p(class = "guide-subtitle", "What is a duty cycle?"),
          tags$ul(
            tags$li("In bioacoustics, a duty cycle is a recording schedule in which the recorder does not run continuously but records at regular intervals — for example 1 min every 5 min, or 10 min every 30 min.")
          ),
          p(class = "guide-subtitle", "Purpose of the tab"),
          tags$ul(
            tags$li("Evaluate how species richness changes under different duty-cycle schedules within the same fixed daily recording window."),
            tags$li("The user defines a range of duty-cycle values to test, from a minimum to a maximum period."),
            tags$li("Compare duty-cycle designs by visualising richness retention and recording effort."),
            tags$li("Identify which species are most likely to be missed under reduced effort.")
          ),
          p(class = "guide-subtitle", "User inputs"),
          tags$ul(
            tags$li(strong("Recording window (left preview): "), "select a start and end event for the recording window (dawn, sunrise, sunset, dusk). Tick '24h recording' if your recorder runs all day."),
            tags$li(strong("Extend before / after: "), "shift the window earlier or later relative to the selected solar events."),
            tags$li(strong("Fixed duration: "), "if used, the window ends after a fixed number of hours rather than at the end event."),
            tags$li(strong("Second daily window: "), "optionally add a second recording window with its own settings."),
            tags$li(strong("Min / Max period (min): "), "the smallest and largest recording interval to test."),
            tags$li(strong("Step (min): "), "the increment between tested period values."),
            tags$li(strong("Duty duration (min): "), "the duration of each individual recording block within a period."),
            tags$li(strong("Bootstrap: "), "number of bootstrap replicates used to assess variability across survey days."),
            tags$li(strong("Show average + CI: "), "overlay the mean richness-effort curve and its 95% bootstrap confidence band."),
            div(style = "color: #A3B3FF; margin-top: 8px; margin-bottom: 8px;",
                "Example: with Min = 1, Max = 30, Step = 5, Duty duration = 1, the app tests: 1 min/1 min, 1 min/6 min, 1 min/11 min … up to 1 min/30 min.")
          ),
          p(class = "guide-subtitle", "Main calculations"),
          tags$ul(
            tags$li("For each duty-cycle schedule, the app simulates subsampling of the detection data."),
            tags$li("At every bootstrap iteration, survey days are resampled with replacement and richness is recalculated."),
            tags$li("Richness is expressed as mean daily richness and standardised to 100% within each replicate."),
            tags$li("Species detected under each schedule are compared with those from the reference (continuous) schedule in the same replicate."),
            tags$li("A species is counted as missed if it appears in the reference but not in the tested schedule.")
          ),
          p(class = "guide-subtitle", "Outputs"),
          tags$ul(
            tags$li(strong("Richness-effort curve: "), "shows how relative richness changes with recording effort. Each curve is one bootstrap replicate; the blue band shows the mean ± 95% CI."),
            tags$li(strong("Missed-species heatmap: "), "shows how often each species goes undetected at each duty-cycle setting compared to continuous recording.")
          ),
          div(
            class = "guide-important",
            p(strong("Important")),
            tags$ul(
              tags$li("A species is counted as missed only if it is detected under the reference schedule but absent under the tested schedule in that same bootstrap replicate.")
            )
          )
        ),
        ###User guide TAB3 ----
        tags$details(
          tags$summary("3. Time window"),
          p(class = "guide-subtitle", "Purpose of the tab"),
          tags$ul(
            tags$li("Evaluate how species richness changes when the recording window is shifted in time and varied in duration."),
            tags$li("Identify which part of the day or night is most efficient for detecting species."),
            tags$li("Compare one-window and two-window designs.")
          ),
          p(class = "guide-subtitle", "User inputs"),
          tags$ul(
            tags$li(strong("Reference event: "), "the anchor for all tested windows (e.g. dawn, sunrise, or a fixed clock time)."),
            tags$li(strong("Min / Max offset: "), "range of start times to test relative to the reference event."),
            tags$li(strong("Start time step: "), "increment between tested start offsets."),
            tags$li(strong("Min / Max duration: "), "range of recording window durations to test."),
            tags$li(strong("Duration step: "), "increment between tested durations."),
            tags$li(strong("Window 2: "), "optionally add a second window with its own reference event, offset range, and duration range."),
            tags$li(
              strong("Keep top combined windows: ")," This parameter limits the cartesian product by retaining,
              for each recording duration, only the N best-performing Window 1 candidates and N best-performing Window 2 candidates before forming combinations.Increase this value if you want a more exhaustive search at the cost of
              longer computation and heavier plots."),
            tags$li(strong("Show optimal frontier: "), "highlight the best-performing design for each duration on the heatmap."),
            tags$li("Note: if Min offset is greater than Max offset the analysis will return an error. Ensure Min offset < Max offset.")
          ),
          p(class = "guide-subtitle", "Main calculations"),
          tags$ul(
            tags$li("A grid of candidate schedules is built from all tested start offsets and durations. The number and parameters of candidate schedules are summarised in the Design preview table."),
            tags$li("For each date, the app simulates subsampling under all candidate schedules and counts distinct species detected."),
            tags$li("If a schedule crosses midnight, detections from the adjacent day are included."),
            tags$li("Mean daily richness is standardised to 100% (= the Window 1 design with the highest mean richness)."),
            tags$li("If both windows are used, the app also evaluates all pairwise combinations. Combined richness is the union of species detected in both windows, and 100% corresponds to the best combined design.")
          ),
          p(class = "guide-subtitle", "Output plots"),
          tags$ul(
            tags$li("Window 1 only: a heatmap of richness across start offsets and durations."),
            tags$li("Window 2 active: a second heatmap for Window 2 and a scatter plot of combined two-window richness versus total effort.")
          ),
          div(
            class = "guide-important",
            p(strong("Important")),
            tags$ul(
              tags$li("This tab identifies efficient recording windows, not duty-cycle designs."),
              tags$li("If tested windows extend beyond the real recording coverage of the dataset, richness may be underestimated."),
              tags$li("When both windows are active, the number of candidate combined designs grows as the product
              of all Window 1 and Window 2 candidates.") ,
              div(style = "color: #A3B3FF; margin-top: 8px; margin-bottom: 8px;","for example 12 start offsets × 12 durations per
                  window produces 144 × 144 = 20 736 combined designs, each of which must be evaluated for
                  every duty-cycle period in Tab 4. A value of 3 reduces the combined set by more
                  than 94% while preserving the designs most likely to appear on the Pareto frontier, since a
                  poorly-performing window is unlikely to produce an optimal combined design regardless of what
                  it is paired with"
              )
            )
          )
        ),
        ###User guide TAB4 ----
        tags$details(
          tags$summary("4. Multi optimum"),
          p(class = "guide-subtitle", "Purpose"),
          tags$ul(
            tags$li("Combine Duty cycle and Time window results into complete sampling designs (daily window + duty cycle)."),
            tags$li("Compare all candidate combinations in terms of richness retention and effort."),
            tags$li("Identify the Pareto-optimal design: highest richness for a given effort.")
          ),
          p(class = "guide-subtitle", "Before running this tab"),
          tags$ul(
            tags$li("The Duty cycle and Time window computations must be run first."),
            tags$li("No additional user input is required; candidate schedules are taken from Tabs 2 and 3.")
          ),
          p(class = "guide-subtitle", "Main calculations"),
          tags$ul(
            tags$li("For each candidate time window × duty-cycle combination, the app builds a final schedule by keeping only the duty-cycle blocks that fall inside the selected window(s)."),
            tags$li("Richness is recalculated on these combined schedules."),
            tags$li("Relative richness is scaled to the maximum richness found across all tested designs.")
          ),
          p(class = "guide-subtitle", "Outputs"),
          tags$ul(
            tags$li(strong("Multi-optimum plot: "), "all tested combinations in the effort–richness space. Recorder battery life equivalents are shown in the tooltip of Pareto-optimal points."),
            tags$li(strong("Pareto frontier: "), "highlighted in red — the set of designs that retain the highest richness for a given recording effort.")
          ),
          div(
            class = "guide-important",
            p(strong("Important")),
            tags$ul(
              tags$li("This tab is the final decision-support step of the application."),
              tags$li("The optimisation criterion is mean daily richness, not cumulative seasonal richness.")
            )
          )
        ),
        ###User guide TAB5 ----
        tags$details(
          tags$summary("5. Target species activity"),
          p(class = "guide-subtitle", "Purpose of the tab"),
          tags$ul(
            tags$li("Visualise the temporal distribution of detections for one or several focal species."),
            tags$li("Assess how different subsampling regimes interact with the timing of species activity.")
          ),
          p(class = "guide-subtitle", "User inputs"),
          tags$ul(
            tags$li("Select one or several target species."),
            tags$li("Choose a range of duty-cycle periods to test (same parameters as Tab 2)."),
            tags$li("Choose the time resolution at which detections are aggregated for display.")
          ),
          p(class = "guide-subtitle", "Main calculations"),
          tags$ul(
            tags$li("For each tested duty-cycle period, the app simulates subsampling of the focal species detections."),
            tags$li("Detections falling into successive time bins across dates are counted.")
          ),
          p(class = "guide-subtitle", "Outputs"),
          tags$ul(
            tags$li(strong("Animated heatmap: "), "shows how apparent activity patterns change under different subsampling regimes.")
          ),
          div(
            class = "guide-important",
            p(strong("Important")),
            tags$ul(
              tags$li("This tab is descriptive rather than optimising."),
              tags$li("It visualises temporal activity patterns but does not identify optimal richness-effort trade-offs.")
            )
          )
        ),
        ##User guide warnings----
        h2("Interpretation cautions"),
        div(
          class = "guide-main",
          tags$ol(
            tags$li(strong("Results depend on the input dataset. "), "The app does not produce universal rules. The optimal design depends on the study site, season, target community, and recording context."),
            tags$li(strong("Richness metrics differ by tab. "), "Most analyses use mean daily richness, not cumulative seasonal richness."),
            tags$li(strong("Missed-species results are strict. "), "A species is counted as missed only if it disappears completely from a bootstrap replicate relative to the reference schedule."),
            tags$li(strong("Small datasets require cautious interpretation. "), "A high bootstrap count does not compensate for a small number of independent recording days.")
          )
        )
      )
    )
  )
)


# SERVER ----

server <- function(input, output, session) {
  
  message(Sys.time(), " session launched")
  rv_filters <- reactiveValues(
    date_active = FALSE,
    date_range  = NULL,
    min_conf    = 0,
    site        = NULL,
    recorder    = NULL,
    year        = NULL,
    species     = NULL
  )
  # Store values between tabs----
  values <- reactiveValues(
    summary_richness    = NULL,
    grid_results        = NULL,
    schedules           = NULL,
    missed_summary      = NULL,
    experimental_domain = NULL
  )
# TAB 1 IMPORT AND METADATA ----
## Read raw BirdNET files  ----
  bn_tables <- eventReactive(input$load_data, {
    
    req(input$bn_file)
    showNotification("Running: loading data...", type = "message", duration = 2)
    
    withProgress(message = "Loading and preparing BirdNET data...", value = 0, {
      
      files <- input$bn_file$datapath
      names <- input$bn_file$name
      
      ### Expand any zip archives ----
      expanded_files <- c()
      for (i in seq_along(files)) {
        if (grepl("\\.zip$", names[i], ignore.case = TRUE)) {
          tmp_dir <- file.path(tempdir(), paste0("unz_", i))
          dir.create(tmp_dir, showWarnings = FALSE)
          unzip(files[i], exdir = tmp_dir)
          zipped_files <- list.files(
            tmp_dir,
            pattern   = "selection\\.table\\.(csv|txt)$",
            full.names = TRUE,
            recursive  = TRUE
          )
          expanded_files <- c(expanded_files, zipped_files)
        } else {
          expanded_files <- c(expanded_files, files[i])
        }
      }
      
      incProgress(0.2, detail = "Reading files...")
      data_list <- lapply(expanded_files, fread)
      
      incProgress(0.6, detail = "Combining files...")
      result <- rbindlist(data_list, fill = TRUE)
      names(result) <- make.names(names(result))
      
      if ("Begin.Path" %in% names(result)) {
        result[, Begin.Path := normalize_begin_path(Begin.Path)]
      }
      
      incProgress(1)
      result
    })
  })
  
  ## Reset path-parser selections when new data is loaded ----
  observeEvent(input$load_data, {
    updateSelectInput(session, "site_col", selected = "None")
    updateSelectInput(session, "rec_col",  selected = "None")
  })
  
  observeEvent(input$apply_date_filter, {
    req(input$date_subset)
    rv_filters$date_active <- TRUE
    rv_filters$date_range  <- input$date_subset
  })
  
  observeEvent(input$reset_date_filter, {
    rv_filters$date_active <- FALSE
    rv_filters$date_range  <- NULL
  })
  ## Extract metadata and standardise column names ----
  dt_parsed <- reactive({
    req(bn_tables())
    
    withProgress(message = "Formatting BirdNET data...", value = 0, {
      
      dt <- bn_tables()
      need_cols(dt, c("Begin.Path", "Common.Name", "Confidence", "Begin.Time..s."),
                context = "BirdNET selection table")
      
      dt[, filename := tools::file_path_sans_ext(basename(Begin.Path))]
      
      incProgress(0.2, detail = "Extracting metadata")
      
      BIRDNET <- dt %>%
        mutate(
          species            = Common.Name,
          conf               = Confidence,
          start              = as.numeric(Begin.Time..s.),
          timestamp          = str_extract(filename, "\\d{8}_\\d{6}"),
          timestamp_fm       = parse_date_time(timestamp, orders = "Ymd_HMS", tz = input$timezone),
          timestamp_adjusted = timestamp_fm + seconds(start)
        ) %>%
        filter(!is.na(timestamp)) %>%
        mutate(
          date           = str_replace(timestamp,
                                       "^(\\d{4})(\\d{2})(\\d{2})_(\\d{2})(\\d{2})(\\d{2})$",
                                       "\\1-\\2-\\3"),
          year           = year(timestamp_adjusted),
          time           = format(timestamp_adjusted, "%H:%M:%S"),
          start_file     = format(timestamp_fm,       "%H:%M:%S"),
          timestamp      = format(timestamp_adjusted, "%Y-%m-%d %H:%M:%S"),
          minute_of_day  = hour(timestamp_adjusted) * 60 + minute(timestamp_adjusted),
          site           = "all",
          recorder       = "global"
        ) %>%
        select(species, site, recorder, filename, start,
               conf, date, year, start_file, timestamp, minute_of_day, Begin.Path)
      
      incProgress(0.8, detail = "Finalizing")
      as.data.table(BIRDNET)
    })
  })
  
  ## Begin.Path preview ----
  output$beginpath_preview <- renderTable({
    req(bn_tables())
    dt <- bn_tables()
    if (!"Begin.Path" %in% names(dt)) return(NULL)
    mat     <- split_path_matrix_preview(dt$Begin.Path)
    preview <- as.data.frame(mat)
    colnames(preview) <- paste0("V", seq_len(ncol(preview)))
    head(preview, 6)
  })
  
  ## path parser electors ----
  output$path_parser_ui <- renderUI({
    req(bn_tables())
    dt <- bn_tables()
    if (!"Begin.Path" %in% names(dt)) return(tags$div("Begin Path column missing"))
    mat <- split_path_matrix_preview(dt$Begin.Path)
    n   <- ncol(mat)
    tagList(
      selectInput("site_col", "Site column",     c("None", paste0("V", seq_len(n)))),
      selectInput("rec_col",  "Recorder column", c("None", paste0("V", seq_len(n))))
    )
  })
  
  ## Apply path parsing  ----
  dt_parsed_path <- reactiveVal(NULL)
  
  observeEvent(input$apply_parsing, {
    req(dt_parsed())
    withProgress(message = "Applying path parsing...", value = 0, {
      
      dt  <- copy(dt_parsed())
      mat <- split_path_matrix(dt$Begin.Path)
      
      site_col <- if (!is.null(input$site_col) && input$site_col != "None")
        as.integer(sub("V", "", input$site_col)) else NULL
      
      rec_col <- if (!is.null(input$rec_col) && input$rec_col != "None")
        as.integer(sub("V", "", input$rec_col)) else NULL
      
      incProgress(0.5)
      dt[, site     := if (!is.null(site_col)) mat[, site_col] else "all"]
      dt[, recorder := if (!is.null(rec_col))  mat[, rec_col]  else "global"]
      
      incProgress(0.5)
      dt_parsed_path(dt)
    })
  })
  
  ## Reset parsed path ----
  observeEvent(input$load_data, { dt_parsed_path(NULL) })
  
  ##parsed path version (if available), else raw parsed ----
  dt_for_app <- reactive({
    req(dt_parsed())
    parsed <- dt_parsed_path()
    if (!is.null(parsed)) return(parsed)
    dt_parsed()
  })
  
  dt_filtered <- reactiveVal(NULL)


  recompute_filters <- function() {
    req(dt_for_app())
    dt <- copy(dt_for_app())
    
    if (isTRUE(rv_filters$date_active) &&
        !is.null(rv_filters$date_range) &&
        length(rv_filters$date_range) == 2) {
      dt <- dt[
        as.Date(date) >= rv_filters$date_range[1] &
          as.Date(date) <= rv_filters$date_range[2]
      ]
    }
    
    if (!is.null(rv_filters$min_conf)) dt <- dt[conf >= rv_filters$min_conf]
    if (!is.null(rv_filters$site))     dt <- dt[site %in% rv_filters$site]
    if (!is.null(rv_filters$recorder)) dt <- dt[recorder %in% rv_filters$recorder]
    if (!is.null(rv_filters$year))     dt <- dt[year %in% rv_filters$year]
    if (!is.null(rv_filters$species))  dt <- dt[species %in% rv_filters$species]
    
    dt_filtered(dt)
  }
  
  
  
# ## Reset filtered dataset when new data is loaded ----
  observeEvent(input$load_data, { dt_filtered(NULL) })

  ## Current working dataset (filtered if applied) ----
  dt_raw <- reactive({
    req(dt_for_app())
    filtered <- dt_filtered()
    if (is.null(filtered)) return(dt_for_app())
    filtered
  })
  
  
  observeEvent(input$apply_date_filter, {
    req(input$date_subset)
    
    withProgress(message = "Filtering dates...", value = 0, {
      
      rv_filters$date_active <- TRUE
      rv_filters$date_range  <- input$date_subset
      
      incProgress(1, detail = "Done")
      recompute_filters()
    })
  })
  observeEvent(input$apply_filters, {
    withProgress(message = "Filtering data...", value = 0, {
      
      rv_filters$min_conf <- input$min_conf
      rv_filters$site     <- input$site
      rv_filters$recorder <- input$recorder
      rv_filters$year     <- input$year
      rv_filters$species  <- input$species
      
      incProgress(1, detail = "Done")
      recompute_filters()
    })
  })
  # observeEvent(input$apply_date_filter, {

  # observeEvent(input$apply_filters, {
  #   req(dt_for_app())
  #   withProgress(message = "Filtering data...", value = 0, {
  #     
  #     dt <- copy(dt_for_app())
  #     
  #     incProgress(0.2, detail = "Filtering by confidence")
  #     if (!is.null(input$min_conf)) dt <- dt[conf >= input$min_conf]
  #     
  #     incProgress(0.2, detail = "Filtering by site / recorder")
  #     if (!is.null(input$site))     dt <- dt[site     %in% input$site]
  #     if (!is.null(input$recorder)) dt <- dt[recorder %in% input$recorder]
  #     
  #     incProgress(0.2, detail = "Filtering by year")
  #     if (!is.null(input$year)) dt <- dt[year %in% input$year]
  #     
  #     incProgress(0.2, detail = "Filtering by species")
  #     if (!is.null(input$species)) dt <- dt[species %in% input$species]
  #     
  #     incProgress(0.2, detail = "Done")
  #     dt_filtered(dt)
  #   })
  # })

 
  ### Filter selector UIs ----
  output$site_ui    <- renderUI({ req(dt_for_app()); selectInput("site",     "Site",    unique(dt_for_app()$site),             multiple = TRUE) })
  output$rec_ui     <- renderUI({ req(dt_for_app()); selectInput("recorder", "Recorder",unique(dt_for_app()$recorder),         multiple = TRUE) })
  output$year_ui    <- renderUI({ req(dt_for_app()); selectInput("year",     "Year",    unique(dt_for_app()$year),             multiple = TRUE) })
  output$species_ui <- renderUI({ req(dt_for_app()); selectInput("species",  "Species", sort(unique(dt_for_app()$species)),    multiple = TRUE) })
  
  observeEvent(input$load_data, {
    if (!is.null(input$site))     updateSelectInput(session, "site",     selected = "all")
    if (!is.null(input$recorder)) updateSelectInput(session, "recorder", selected = "global")
  })
  
  ## Data preview table ----
  output$table_preview <- renderTable({
    req(dt_raw())
    head(dt_raw(), 20)
  })
  
  
  ## Coordinates ----
  ### multiple coordinates ----
  coord_multiple <- reactive({
    req(input$coord_mode == "multiple", input$coord_file)
    read_coord_file(input$coord_file$datapath, input$coord_file$name)
  })
  
  ### Single coordinates ----
  active_lat <- reactive({
    if (input$coord_mode == "single") input$single_lat
    else { req(coord_multiple()); coord_multiple()$latitude[1] }
  })
  
  active_lon <- reactive({
    if (input$coord_mode == "single") input$single_lon
    else { req(coord_multiple()); coord_multiple()$longitude[1] }
  })
  
  ##Solar event ----
  ### Solar for single pair coord ----
  solar_table <- reactive({
    req(dt_for_app(), input$timezone)
    dt          <- dt_for_app()
    dates_unique <- unique(dt$date)
    
    if (input$coord_mode == "single") {
      # One coordinate pair for all dates ----
      lat <- input$single_lat
      lon <- input$single_lon
      solar_list <- lapply(dates_unique, function(d) {
        sun <- getSunlightTimes(date = as.Date(d), lat = lat, lon = lon,
                                keep = c("dawn","sunrise","sunset","dusk"), tz = input$timezone)
        data.table(date = as.Date(d), recorder = "global",
                   dawn = sun$dawn, sunrise = sun$sunrise, sunset = sun$sunset, dusk = sun$dusk)
      })
      return(rbindlist(solar_list))
    }
    
    ### Solar for multiple coord ----
    req(coord_multiple())
    coords <- coord_multiple()
    req("recorder" %in% names(coords), "latitude" %in% names(coords), "longitude" %in% names(coords))
    
    comb <- unique(dt[, .(recorder, date)])
    comb <- merge(comb, coords, by = "recorder", all.x = TRUE)
    
    solar_list <- lapply(seq_len(nrow(comb)), function(i) {
      sun <- getSunlightTimes(date = as.Date(comb$date[i]),
                              lat  = comb$latitude[i],
                              lon  = comb$longitude[i],
                              keep = c("dawn","sunrise","sunset","dusk"), tz = input$timezone)
      data.table(recorder = comb$recorder[i], date = as.Date(comb$date[i]),
                 dawn = sun$dawn, sunrise = sun$sunrise, sunset = sun$sunset, dusk = sun$dusk)
    })
    rbindlist(solar_list)
  })
  
  # TAB 2 DUTY CYCLE ----
  ## Schedule preview plots ----
  
  ### Left panel: daily recording window ----
  output$sun_plot <- renderPlot({
    req(input$start_event, input$end_event, input$extend_before, input$extend_after)
    showNotification("Schedule preview updated", type = "message", duration = 1)
    
    sun_events <- getSunlightTimes(
      date = Sys.Date(), lat = active_lat(), lon = active_lon(),
      keep = c("dawn","sunrise","sunset","dusk"), tz = input$timezone
    )
    sun_events <- as.data.table(sun_events)
    sun_events <- melt(sun_events, measure.vars = c("dawn","sunrise","sunset","dusk"),
                       variable.name = "event", value.name = "datetime")
    setDT(sun_events)
    sun_events[, minute_of_day := hour(datetime) * 60 + minute(datetime)]
    
    #### Window 1 boundaries ----
    start_min1 <- sun_events[event == input$start_event, minute_of_day]
    if (input$length_morning > 0) {
      end_min1 <- start_min1 + input$length_morning * 60
    } else {
      end_min1 <- sun_events[event == input$end_event, minute_of_day]
    }
    start_min1 <- (start_min1 - input$extend_before * 60) %% (24 * 60)
    end_min1   <- (end_min1   + input$extend_after  * 60) %% (24 * 60)
    
    if (isTRUE(input$record_24h)) { start_min1 <- 0; end_min1 <- 1440 }
    
    #### Window 2 boundaries ----

    if (isTRUE(input$use_window2)) {
      start_min2 <- sun_events[event == input$start_event2, minute_of_day]
      if (input$length_morning2 > 0) {
        end_min2 <- start_min2 + input$length_morning2 * 60
      } else {
        end_min2 <- sun_events[event == input$end_event2, minute_of_day]
      }
      start_min2 <- (start_min2 - input$extend_before2 * 60) %% (24 * 60)
      end_min2   <- (end_min2   + input$extend_after2  * 60) %% (24 * 60)
    }
    
    #### Helper to add a window rectangle----
    add_window_rect <- function(p, s, e, fill_col) {
      if (s <= e) {
        p + geom_rect(aes(xmin = s, xmax = e), ymin = 0.9, ymax = 1.1, fill = fill_col, alpha = 0.5)
      } else {
        p +
          geom_rect(aes(xmin = s,    xmax = 1440), ymin = 0.9, ymax = 1.1, fill = fill_col, alpha = 0.5) +
          geom_rect(aes(xmin = 0,    xmax = e),    ymin = 0.9, ymax = 1.1, fill = fill_col, alpha = 0.5)
      }
    }
    
    ###Left panel plot----
    p <- ggplot()
    p <- add_window_rect(p, start_min1, end_min1, "skyblue")
    if (isTRUE(input$use_window2)) p <- add_window_rect(p, start_min2, end_min2, "orange")
    
    p +
      geom_vline(data = sun_events,
                 aes(xintercept = minute_of_day, color = event), linewidth = 0.8) +
      scale_color_manual(values = c(dawn = "darkorange", sunrise = "gold",
                                    sunset = "firebrick", dusk = "purple")) +
      scale_x_continuous(limits = c(0, 1440), breaks = NULL) +
      theme_minimal() +
      theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), legend.position = "bottom", legend.title = element_blank())
  })
  
  ### Right panel: duty-cycle blocks----
  output$duty_plot <- renderPlot({
    req(input$period_min, input$period_max, input$period_step, input$duty_duration)
    
    periods       <- seq(input$period_min, input$period_max, by = input$period_step)
    duty_duration <- input$duty_duration
    
    all_duty_cycles <- rbindlist(
      lapply(periods, function(p) {
        start_times <- seq(0, 1440 - duty_duration, by = p)
        data.table(period = p, start = start_times, end = start_times + duty_duration)
      })
    )
    
    # Show only the morning half (0–720 min) for readability 
    dt_am <- all_duty_cycles[start < 720 & end > 0]
    dt_am[, start_plot := pmax(start, 0)]
    dt_am[, end_plot   := pmin(end,   720)]
    dt_am <- dt_am[start_plot < end_plot]
    
    #### Right panel plot----
    base_theme <- theme_minimal(base_size = 11) +
      theme(axis.text.y = element_text(size = 9, face = "bold"),
            axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
            panel.grid.major.x = element_line(color = "grey90", linewidth = 0.5),
            panel.grid.minor.x = element_line(color = "grey95", linewidth = 0.3),
            panel.grid.major.y = element_blank(), axis.title.x = element_blank(),
            plot.margin = margin(10, 10, 10, 10))
    
    ggplot(dt_am) +
      geom_segment(
        aes(x = start_plot, xend = end_plot,
            y    = factor(period, levels = rev(unique(period))),
            yend = factor(period, levels = rev(unique(period)))),
        color = "skyblue", alpha = 0.8, linewidth = 4.5, lineend = "butt"
      ) +
      scale_x_continuous(limits = c(0, 720), breaks = seq(0, 720, by = 60),
                         labels = function(x) sprintf("%02d:00", x / 60)) +
      scale_y_discrete(name = "Period (min)", labels = function(x) paste(x, "min")) +
      base_theme
  })
  
  
  # TAB2: Duty-cycle  ----
  
  temporal_results <- eventReactive(input$run_temporal, {
    
    req(input$period_step > 0)
    message(Sys.time(), " temporal rarefaction: start")
    
    periods    <- seq(as.numeric(input$period_min), as.numeric(input$period_max),
                      by = as.numeric(input$period_step))
    param_grid <- CJ(period = periods, duty_duration = input$duty_duration)
    
    local_dt <- dt_raw()
    req(local_dt)
    
    ## Build schedules ----
    withProgress(message = "Creating schedules...", value = 0, {
      
      recs     <- unique(local_dt$recorder)
      solar_dt <- setDT(solar_table())
      
      schedules <- lapply(seq_len(nrow(param_grid)), function(i) {
        
        p <- param_grid[i]
        
        sched_rec_list <- lapply(recs, function(r) {
          sun_r <- solar_dt[recorder == r]
          
          sched_r1 <- create_sun_duty_schedule(
            sun_table           = sun_r,
            start_event         = input$start_event,
            end_event           = input$end_event,
            length_morning      = input$length_morning,
            extend_before_hours = input$extend_before,
            extend_after_hours  = input$extend_after,
            period              = p$period,
            duty_duration       = p$duty_duration,
            nb_duty             = 1,
            random_start        = FALSE,
            record_24h          = input$record_24h
          )
          
          if (isTRUE(input$use_window2)) {
            sched_r2 <- create_sun_duty_schedule(
              sun_table           = sun_r,
              start_event         = input$start_event2,
              end_event           = input$end_event2,
              length_morning      = input$length_morning2,
              extend_before_hours = input$extend_before2,
              extend_after_hours  = input$extend_after2,
              period              = p$period,
              duty_duration       = p$duty_duration,
              nb_duty             = 1,
              random_start        = FALSE
            )
            sched_r <- rbindlist(list(sched_r1, sched_r2))
          } else {
            sched_r <- sched_r1
          }
          
          sched_r[, recorder := r]
          sched_r
        })
        
        sched <- rbindlist(sched_rec_list)
        setattr(sched, "period",        p$period)
        setattr(sched, "duty_duration", p$duty_duration)
        sched
      })
      
      values$schedules <- schedules
    })
    
    ## Store experimental design metadata for Tab 4 ----
    values$experimental_domain <- list(
      timezone = input$timezone,
      w1 = list(start_event = input$start_event, end_event = input$end_event,
                length = input$length_morning, extend_before = input$extend_before,
                extend_after = input$extend_after),
      w2 = if (isTRUE(input$use_window2)) list(
        start_event = input$start_event2, end_event = input$end_event2,
        length = input$length_morning2, extend_before = input$extend_before2,
        extend_after = input$extend_after2) else NULL
    )
    
    ## Bootstrap loop ----
    withProgress(message = "Running bootstrap iterations...", value = 0, {
      
      all_results <- vector("list", input$B_temporal)
      
      if (!is.data.table(local_dt)) local_dt <- as.data.table(local_dt)
      local_dt[, date := as.Date(date)]
      all_dates  <- unique(local_dt$date)
      n_days     <- length(all_dates)
      
      if (length(values$schedules) == 0) stop("No schedules found.")
      
      schedule_periods <- sapply(values$schedules, function(s) attr(s, "period"))
      n_schedules      <- length(values$schedules)
      
      ## sum total effort ----
      filtered_dates_dt <- data.table(date = unique(as.Date(local_dt$date)))
      schedule_efforts  <- sapply(values$schedules, function(s) {
        s_sub <- merge(as.data.table(s), filtered_dates_dt, by = "date", all = FALSE)
        sum(s_sub$end_min - s_sub$start_min + 1)
      })
      
      for (b in seq_len(input$B_temporal)) {
        
        incProgress(1 / input$B_temporal, detail = paste(b, "/", input$B_temporal))
        
        sampled_days   <- sample(all_dates, n_days, replace = TRUE)
        sampled_days_dt <- data.table(date = as.Date(sampled_days))
        dt_boot <- local_dt[sampled_days_dt, on = "date", allow.cartesian = TRUE, nomatch = 0L]
        
        if (nrow(dt_boot) == 0) {
          all_results[[b]] <- data.table(
            period           = schedule_periods,
            effort           = schedule_efforts,
            richness         = 0,
            bootstrap        = b,
            species_detected = list(rep(list(character(0)), n_schedules))
          )
          next
        }
        
        one_bootstrap <- vector("list", n_schedules)
        
        for (i in seq_len(n_schedules)) {
          
          s   <- values$schedules[[i]]
          res <- tryCatch(
            compute_richness(dt = dt_boot, schedule = s,
                             per_spatial = "all", per_temporal = "day",
                             return_species = TRUE),
            error = function(e) NULL
          )
          
          r_mean      <- 0
          sp_detected <- character(0)
          
          if (!is.null(res) && !is.null(res$richness) && nrow(res$richness) > 0) {
            
            res_dt   <- as.data.table(res$richness)
            date_col <- if ("date" %in% names(res_dt)) "date" else
              names(res_dt)[which(sapply(res_dt, inherits, "Date"))[1]]
            val_col  <- if ("richness" %in% names(res_dt)) "richness" else
              names(res_dt)[which(sapply(res_dt, is.numeric) & names(res_dt) != date_col)[1]]
            
            if (!is.null(date_col) && !is.null(val_col)) {
              r_subset <- res_dt[, list(date = get(date_col), richness = get(val_col))]
              r_merged <- merge(sampled_days_dt, r_subset, by = "date", all.x = TRUE)
              r_merged[is.na(richness), richness := 0]
              r_mean   <- mean(r_merged$richness)
              
              if (!is.null(res$species) && nrow(res$species) > 0) {
                sp_dt   <- as.data.table(res$species)
                sp_col  <- if ("species" %in% names(sp_dt)) "species" else
                  names(sp_dt)[which(sapply(sp_dt, function(x) is.character(x) || is.factor(x)))[1]]
                if (!is.null(sp_col)) {
                  if ("date" %in% names(sp_dt))
                    sp_detected <- unique(sp_dt[date %in% sampled_days, get(sp_col)])
                  else
                    sp_detected <- unique(sp_dt[[sp_col]])
                }
              }
            }
          }
          
          one_bootstrap[[i]] <- data.table(
            period           = schedule_periods[i],
            effort           = schedule_efforts[i],
            richness         = r_mean,
            bootstrap        = b,
            species_detected = list(sp_detected)
          )
        }
        
        all_results[[b]] <- rbindlist(one_bootstrap)
      }
      
      summary_richness <- rbindlist(all_results)
      summary_richness[, richness_pct        := 100 * richness / max(richness), by = bootstrap]
      summary_richness[, n_species_detected  := lengths(species_detected)]
      
    ## Missed-species analysis ----
      min_period_val  <- min(summary_richness$period)
      reference_table <- summary_richness[
        period == min_period_val,
        list(reference_species = list(unique(unlist(species_detected)))),
        by = "bootstrap"
      ]
      
      summary_richness <- merge(summary_richness, reference_table, by = "bootstrap", all.x = TRUE)
      
      missed_species_list <- summary_richness[, list(
        species_missed = mapply(
          function(ref, detected) setdiff(ref, unlist(detected)),
          reference_species, species_detected, SIMPLIFY = FALSE
        )
      ), by = list(bootstrap, period)]
      
      missed_species_long <- missed_species_list[,
                                                 list(species_missed = unlist(species_missed)),
                                                 by = list(bootstrap, period)]
      missed_summary <- missed_species_long[, .N, by = list(species_missed, period)]
      
      species_order <- missed_summary[
        , list(total_missed = sum(N)), by = "species_missed"
      ][order(total_missed), species_missed]
      
      missed_summary[, species_missed := factor(species_missed, levels = species_order)]
      
      # Fill all (species, period) combinations so heatmap has no gaps
      all_periods <- sort(unique(summary_richness$period))
      missed_summary_full <- CJ(species_missed = species_order, period = all_periods)
      missed_summary_full <- merge(missed_summary_full, missed_summary,
                                   by = c("species_missed", "period"), all.x = TRUE)
      missed_summary_full[is.na(N), N := 0]
      missed_summary_full[, species_missed := factor(species_missed, levels = species_order)]
      
      values$missed_summary    <- missed_summary_full
      values$summary_richness  <- summary_richness
      
      message(Sys.time(), " temporal rarefaction: done")
      summary_richness
    })
  })
  
  observeEvent(input$return_duty_data, {
    assign("duty_cycle_data", temporal_results(), envir = .GlobalEnv)
  })
  
  
  # TAB 2 PLOTS ----
  ##bootstrap plot ----
  gg_raref <- reactive({
    req(temporal_results())
    validate(
      need(nrow(temporal_results()) > 0,
           "Invalid temporal window.\nEnd time occurs before start time.\nCheck start/end events or extensions.")
    )
    
    ci_summary <- temporal_results()[, .(
      mean_richness = mean(richness_pct, na.rm = TRUE),
      lower         = quantile(richness_pct, 0.025, na.rm = TRUE),
      upper         = quantile(richness_pct, 0.975, na.rm = TRUE)
    ), by = effort]
    
    gg <- ggplot(
      temporal_results(),
      aes(x = effort, y = richness_pct, group = bootstrap, color = factor(bootstrap))
    ) +
      geom_line(alpha = 0.5, aes(text = paste0(
        "Effort: ", .data$effort, " min — ", round(.data$effort / 60, 2), " h",
        "<br>Richness: ", round(.data$richness_pct, 1), "%",
        "<br>Species detected: ", .data$n_species_detected,
        "<br>Bootstrap ", bootstrap))) +
      geom_point(alpha = 0.5, aes(text = paste0(
        "Effort: ", .data$effort, " min — ", round(.data$effort / 60, 2), " h",
        "<br>Richness: ", round(.data$richness_pct, 1), "%",
        "<br>Species detected: ", .data$n_species_detected,
        "<br>Bootstrap ", bootstrap))) +
      labs(x = "Sampling schedule effort (min)", y = "Relative richness [%]") +
      theme_minimal() +
      guides(color = "none")
    
    if (input$show_summary) {
      gg <- gg +
        geom_ribbon(data = ci_summary,
                    aes(x = .data$effort, ymin = .data$lower, ymax = .data$upper),
                    inherit.aes = FALSE, fill = "darkblue", alpha = 0.4) +
        geom_line(data = ci_summary,
                  aes(x = .data$effort, y = .data$mean_richness),
                  inherit.aes = FALSE, color = "red")
    }
    gg
  })
  
  plotly_raref <- reactive({ ggplotly(gg_raref(), tooltip = "text") })
  output$plot_temporal <- renderPlotly({ plotly_raref() })
  
  observeEvent(input$return_raref, {
    assign("ggplot_rarefaction", gg_raref(),    envir = .GlobalEnv)
    assign("plotly_rarefaction", plotly_raref(), envir = .GlobalEnv)
  })
  
  
  ##Missed-species heatmap ----
  
  gg_missed <- reactive({
    req(values$missed_summary)
    missed_summary <- values$missed_summary
    period_min  <- min(missed_summary$period)
    period_max  <- max(missed_summary$period)
    period_step <- unique(diff(sort(unique(missed_summary$period))))[1]
    
    ggplot(missed_summary, aes(x = period, y = species_missed, fill = N)) +
      geom_tile(aes(text = paste0(
        "Species: ", species_missed, "\nPeriod: ", period, " min\nTimes missed: ", N
      )), color = "white", width = period_step, height = 1) +
      scale_fill_viridis_c(option = "D", name = "# bootstraps missed",
                           limits = c(0, max(missed_summary$N, na.rm = TRUE))) +
      scale_x_continuous(name = "Period (min)",
                         breaks = seq(period_min, period_max, by = period_step),
                         expand = c(0, 0)) +
      labs(y = "Species", title = "Species missed across duty cycles") +
      theme_minimal(base_size = 12) +
      theme(axis.text.y = element_text(size = 8),
            axis.text.x = element_text(angle = 90, vjust = 0.5),
            panel.grid   = element_blank()) +
      geom_vline(xintercept = seq(period_min - period_step / 2,
                                  period_max + period_step / 2,
                                  by = period_step),
                 color = "grey90", linewidth = 0.2) +
      geom_hline(yintercept = seq(0.5, length(levels(missed_summary$species_missed)) + 0.5, by = 1),
                 color = "grey90", linewidth = 0.2)
  })
  
  plotly_missed <- reactive({
    ggplotly(gg_missed(), tooltip = "text") %>%
      plotly::config(modeBarButtonsToRemove = c("select2d","lasso2d","zoomIn2d","zoomOut2d"),
                     displaylogo = FALSE)
  })
  
  output$heatmap_missed_species <- renderPlotly({ plotly_missed() })
  
  observeEvent(input$return_heatmap_miss, {
    assign("gg_heatmap_missed",    gg_missed(),    envir = .GlobalEnv)
    assign("plotly_heatmap_missed", plotly_missed(), envir = .GlobalEnv)
  })
  
  # TAB 3 TIME WINDOW ----
 ## Design preview table ----
  
  output$grid_preview <- DT::renderDataTable({
    
    req(input$grid_start_min, input$grid_start_max, input$grid_start_step,
        input$grid_dur_min,   input$grid_dur_max,   input$grid_dur_step)
    
    starts <- seq(input$grid_start_min, input$grid_start_max, by = input$grid_start_step)
    durs   <- seq(input$grid_dur_min,   input$grid_dur_max,   by = input$grid_dur_step)
    
    design <- as.data.table(expand.grid(start_offset_h = starts, duration_h = durs))
    design[, end_offset_h := start_offset_h + duration_h]
    design[, Window  := sprintf("%+.1fh \u2192 %+.1fh", start_offset_h, end_offset_h)]
    design[, Effort  := sprintf("%.1f h", duration_h)]
    design[, Start   := sprintf("Start = Reference event %+.1f h", start_offset_h)]
    
    if (isTRUE(input$use_window2_grid)) {
      req(input$grid_start_min2, input$grid_start_max2, input$grid_start_step2,
          input$grid_dur_min2,   input$grid_dur_max2,   input$grid_dur_step2)
      starts2 <- seq(input$grid_start_min2, input$grid_start_max2, by = input$grid_start_step2)
      durs2   <- seq(input$grid_dur_min2,   input$grid_dur_max2,   by = input$grid_dur_step2)
      design2 <- as.data.table(expand.grid(start_offset_h2 = starts2, duration_h2 = durs2))
      design2[, end_offset_h2 := start_offset_h2 + duration_h2]
      design2[, Window2 := sprintf("%+.1fh \u2192 %+.1fh", start_offset_h2, end_offset_h2)]
      design2[, Effort2 := sprintf("%.1f h", duration_h2)]
      design <- cbind(design, design2[seq_len(nrow(design))])
    }
    
    setorder(design, start_offset_h, duration_h)
    
    out <- if (isTRUE(input$use_window2_grid))
      design[, .(start_group = Start, Window, Effort, Window2, Effort2)]
    else
      design[, .(start_group = Start, Window, Effort)]
    
    DT::datatable(
      out, rownames = FALSE, extensions = "RowGroup",
      options = list(
        paging = FALSE, scrollY = "150px", scrollCollapse = TRUE,
        info = FALSE, searching = FALSE, ordering = FALSE, dom = "t",
        rowGroup = list(dataSrc = 0),
        columnDefs = list(list(visible = FALSE, targets = 0),
                          list(className = "dt-center", targets = "_all")),
        class = "compact cell-border"
      )
    ) %>% formatStyle(names(out), color = "white")
  })
  
  output$grid_summary <- renderUI({
    req(input$grid_start_min, input$grid_start_max, input$grid_start_step,
        input$grid_dur_min,   input$grid_dur_max,   input$grid_dur_step)
    
    n_starts <- length(seq(input$grid_start_min, input$grid_start_max, by = input$grid_start_step))
    n_durs   <- length(seq(input$grid_dur_min,   input$grid_dur_max,   by = input$grid_dur_step))
    n_total  <- n_starts * n_durs
    
    if (isTRUE(input$use_window2_grid)) {
      n_starts2 <- length(seq(input$grid_start_min2, input$grid_start_max2, by = input$grid_start_step2))
      n_durs2   <- length(seq(input$grid_dur_min2,   input$grid_dur_max2,   by = input$grid_dur_step2))
      n_total   <- n_total + n_starts2 * n_durs2
    }
    
    HTML(sprintf(
      "<small><b>%d designs generated</b> (%d start offsets \u00d7 %d durations) \u2014 reference: <b>%s</b></small>",
      n_total, n_starts, n_durs, input$grid_event
    ))
  })
  
  
  ## Time-window computation ----
  
  grid_results <- eventReactive(input$run_grid, {
    
    dt <- copy(dt_raw())
    req(dt)
    
    withProgress(message = "Building grid...", value = 0, {
      
      dt[, minute := hour(timestamp) * 60 + minute(timestamp)]
      dt[, date   := as.Date(date)]
      setkey(dt, recorder, date)
      
      comb      <- unique(dt[, .(recorder, date)])
      solar_dt  <- setDT(solar_table())
      solar_sub <- solar_dt[, .(recorder, date, dawn, sunrise, sunset, dusk)]
      comb_solar <- solar_sub[comb, on = c("recorder", "date")]
      
      # Convert start-offset and duration inputs from hours to minutes
      start_offsets_min <- seq(input$grid_start_min * 60,
                               input$grid_start_max * 60,
                               by = input$grid_start_step * 60)
      durations_min     <- seq(input$grid_dur_min   * 60,
                               input$grid_dur_max   * 60,
                               by = input$grid_dur_step * 60)
      
      ### Reference minute (W1) ----
      if (input$grid_event %in% c("dawn","sunrise","sunset","dusk")) {
        comb_solar[, ref_min := as.integer(hour(get(input$grid_event)) * 60 +
                                             minute(get(input$grid_event)))]
      } else {
        h_part <- as.integer(sub(":.*", "", input$grid_event))
        m_part <- as.integer(sub(".*:", "", input$grid_event))
        comb_solar[, ref_min := h_part * 60L + m_part]
      }
      
      ## Reference minute (W2) ----
      if (isTRUE(input$use_window2_grid)) {
        if (input$grid_event2 %in% c("dawn","sunrise","sunset","dusk")) {
          comb_solar[, ref_min2 := as.integer(hour(get(input$grid_event2)) * 60 +
                                                minute(get(input$grid_event2)))]
        } else {
          h_part2 <- as.integer(sub(":.*", "", input$grid_event2))
          m_part2 <- as.integer(sub(".*:", "", input$grid_event2))
          comb_solar[, ref_min2 := h_part2 * 60L + m_part2]
        }
      }
      
      ### W1: build grid  ----
      grid_all <- comb_solar[, {
        as.data.table(expand.grid(start_offset_min = start_offsets_min,
                                  duration         = durations_min))
      }, by = .(recorder, date, ref_min)]
      
      grid_all[, `:=`(start_abs_min = ref_min + start_offset_min, effort = duration)]
      
      # Helper: extract detected species for one window 
      species_in_window <- function(rec_i, d_i, s0_i, s1_i) {
        if (s0_i >= 0 && s1_i <= 1440) {
          unique(dt[.(rec_i, d_i)][minute >= s0_i & minute <= s1_i, species])
        } else if (s1_i > 1440) {
          sp1 <- unique(dt[.(rec_i, d_i)][minute >= max(0L, s0_i) & minute <= 1440L, species])
          sp2 <- unique(dt[.(rec_i, as.Date(d_i) + 1L)][minute >= 0L & minute <= (s1_i - 1440L), species])
          unique(c(sp1, sp2))
        } else if (s0_i < 0) {
          sp0 <- unique(dt[.(rec_i, as.Date(d_i) - 1L)][minute >= (1440L + s0_i) & minute <= 1440L, species])
          sp1 <- unique(dt[.(rec_i, d_i)][minute >= 0L & minute <= min(1440L, s1_i), species])
          unique(c(sp0, sp1))
        } else character(0)
      }
      
      grid_all[, species_w1 := mapply(
        species_in_window,
        recorder, date, start_abs_min, start_abs_min + duration,
        SIMPLIFY = FALSE
      )]
      grid_all[, richness := sapply(species_w1, length)]
      
      incProgress(0.4, detail = "Summarising W1...")
      
     
      grid_all[, start_off_h  := start_offset_min / 60]
      grid_all[, duration_h   := duration / 60]
      
    
      grid_daily <- grid_all[, .(
        richness = length(unique(unlist(species_w1))),
        effort   = mean(effort, na.rm = TRUE),
        species_union = list(unique(unlist(species_w1)))   # on garde la liste pour l'étape 2
      ), by = .(date, start_off_h, duration_h)]
      
      # Étape 2 : moyenne sur les jours + union globale des espèces détectées
      grid_w1 <- grid_daily[, .(
        richness           = mean(richness, na.rm = TRUE),
        effort             = mean(effort,   na.rm = TRUE),
        n_species_detected = length(unique(unlist(species_union)))
      ), by = .(start_off_h, duration_h)]
      
      grid_w2       <- NULL
      grid_combined <- NULL
      
      ### W2 build grid----
      if (isTRUE(input$use_window2_grid)) {
        
        start_offsets_min2 <- seq(input$grid_start_min2 * 60,
                                  input$grid_start_max2 * 60,
                                  by = input$grid_start_step2 * 60)
        durations_min2     <- seq(input$grid_dur_min2   * 60,
                                  input$grid_dur_max2   * 60,
                                  by = input$grid_dur_step2 * 60)
        
        grid_all2 <- comb_solar[, {
          as.data.table(expand.grid(start_offset_min = start_offsets_min2,
                                    duration         = durations_min2))
        }, by = .(recorder, date, ref_min = ref_min2)]
        
        grid_all2[, `:=`(start_abs_min = ref_min + start_offset_min, effort = duration)]
        
        grid_all2[, species_w2 := mapply(
          species_in_window,
          recorder, date, start_abs_min, start_abs_min + duration,
          SIMPLIFY = FALSE
        )]
        grid_all2[, richness   := sapply(species_w2, length)]
        grid_all2[, start_off_h := start_offset_min / 60]
        grid_all2[, duration_h  := duration / 60]
        
        grid_w2 <- grid_all2[, .(
          richness           = mean(richness, na.rm = TRUE),
          effort             = mean(effort,   na.rm = TRUE),
          n_species_detected = length(unique(unlist(species_w2)))
        ), by = .(start_off_h, duration_h)]
        
        incProgress(0.7, detail = "Summarising W2...")
        
        ### Combined W1 × W2 ----
        # TOP N windows combined to keep! can me m0dified!
        top_n <- as.integer(input$top_n_combined)
        
        best_w1_keys <- grid_all[, .(
          richness = length(unique(unlist(species_w1)))
        ), by = .(start_off_h, duration_h)][
          order(duration_h, -richness)
        ][, .SD[seq_len(min(.N, top_n))], by = duration_h]
        
        best_w2_keys <- grid_all2[, .(
          richness = length(unique(unlist(species_w2)))
        ), by = .(start_off_h, duration_h)][
          order(duration_h, -richness)
        ][, .SD[seq_len(min(.N, top_n))], by = duration_h]
        
        grid_all_sp  <- grid_all[best_w1_keys, on = c("start_off_h","duration_h"), nomatch = 0][
          , .(recorder, date, start_off_h, duration_h, dur_w1 = duration, species_w1)]
        
        grid_all2_sp <- grid_all2[best_w2_keys, on = c("start_off_h","duration_h"), nomatch = 0][
          , .(recorder, date, start_off_h2 = start_off_h, duration_h2 = duration_h, dur_w2 = duration, species_w2)]
        
        grid_all_comb <- merge(grid_all_sp, grid_all2_sp, by = c("recorder","date"), allow.cartesian = TRUE)
      
        grid_all_comb[, `:=`(
          richness_w1 = sapply(species_w1, function(x) length(unique(unlist(x)))),
          richness_w2 = sapply(species_w2, function(x) length(unique(unlist(x)))),
          richness    = mapply(function(sp1, sp2) length(unique(c(unlist(sp1), unlist(sp2)))),
                               species_w1, species_w2, SIMPLIFY = TRUE),
          effort      = dur_w1 + dur_w2
        )]
        
        grid_combined <- grid_all_comb[, .(
          richness            = mean(richness,    na.rm = TRUE),
          richness_w1         = mean(richness_w1, na.rm = TRUE),
          richness_w2         = mean(richness_w2, na.rm = TRUE),
          effort              = mean(effort,       na.rm = TRUE),
          n_species_detected  = length(unique(c(unlist(species_w1), unlist(species_w2))))
        ), by = .(start_off_h, duration_h, start_off_h2, duration_h2)]
        
        # Clean up large list columns ----
        grid_all_comb[, c("species_w1","species_w2") := NULL]
        
        ### Common normalisation W1, W2, Combined ----
        max_global <- max(c(grid_w1$richness, grid_w2$richness, grid_combined$richness), na.rm = TRUE)
        grid_w1[,       richness_pct    := 100 * richness    / max_global]
        grid_w2[,       richness_pct    := 100 * richness    / max_global]
        grid_combined[, richness_pct    := 100 * richness    / max_global]
        grid_combined[, richness_w1_pct := 100 * richness_w1 / max_global]
        grid_combined[, richness_w2_pct := 100 * richness_w2 / max_global]
        
      } else {
        #### W1 only ----
        grid_w1[, richness_pct := 100 * richness / max(grid_w1$richness, na.rm = TRUE)]
      }
      
      incProgress(1, detail = "Done")
      message(Sys.time(), " window tab: results finalised")
      
      values$grid_results <- list(w1 = grid_w1, w2 = grid_w2, comb = grid_combined)
      values$grid_results
    })
  })
  
  
  # TAB 3 PLOTS ----
  ## heatmap W1 ----
  gg_plot_w1 <- reactive({
    res <- grid_results()
    req(res$w1)
    df <- copy(res$w1)
    
    frontier <- df[, .SD[which.max(richness_pct)], by = duration_h]
    
    dx_vals <- diff(sort(unique(df$start_off_h)))
    dy_vals <- diff(sort(unique(df$duration_h)))
    dx <- if (length(dx_vals) > 0) min(dx_vals) else 0.5
    dy <- if (length(dy_vals) > 0) min(dy_vals) else 0.5
    
    xmin <- min(df$start_off_h) - dx/2;  xmax <- max(df$start_off_h) + dx/2
    ymin <- min(df$duration_h)  - dy/2;  ymax <- max(df$duration_h)  + dy/2
    grid_v <- data.frame(x = unique(df$start_off_h) - dx/2)
    grid_h <- data.frame(y = unique(df$duration_h)  - dy/2)
    
    df[, tooltip := paste0(
      "<b>Window 1</b>",
      "<br>Start offset: ", round(start_off_h, 2), " h",
      "<br>Duration: ",     round(duration_h,  2), " h",
      "<br>Richness: ",     round(richness_pct, 1), " %",
      "<br>Species: ",      n_species_detected)]
    
    p <- ggplot(df, aes(x = start_off_h, y = duration_h)) +
      geom_contour(aes(z = richness_pct), binwidth = 5, color = "grey89", linewidth = 0.2) +
      geom_tile(aes(fill = richness_pct, text = tooltip), color = "grey89", linewidth = 0.5) +
      scale_fill_scico(palette = "oslo", name = "Mean richness (%)", direction = -1) +
      labs(x = paste0("Start offset relative to ", input$grid_event, " (hours)"),
           y = "Recording duration (hours)", title = "Window 1") +
      theme_minimal() +
      theme(panel.grid = element_blank(), axis.ticks = element_line(),
            axis.ticks.length = unit(3, "pt")) +
      scale_x_continuous(breaks = seq(floor(min(df$start_off_h)),   ceiling(max(df$start_off_h)),   by = 0.5)) +
      scale_y_continuous(breaks = seq(floor(min(df$duration_h)),    ceiling(max(df$duration_h)),    by = 0.5)) +
      geom_segment(data = grid_v, aes(x = x, xend = x, y = ymin, yend = ymax),
                   inherit.aes = FALSE, color = "grey89", linewidth = 0.1) +
      geom_segment(data = grid_h, aes(x = xmin, xend = xmax, y = y, yend = y),
                   inherit.aes = FALSE, color = "grey89", linewidth = 0.1)
    
    if (isTRUE(input$show_frontier)) {
      p <- p + geom_point(data = frontier,
                          aes(start_off_h, duration_h,
                              text = paste0("<b>Optimal frontier</b><br>",
                                            "Start offset: ", round(start_off_h, 2), " h<br>",
                                            "Duration: ", round(duration_h, 2), " h<br>",
                                            "Richness: ", round(richness_pct, 1), " %")),
                          inherit.aes = FALSE, color = "brown1", size = 3, shape = 8)
    }
    p
  })
  
  plotly_w1 <- reactive({ ggplotly(gg_plot_w1(), tooltip = "text") })
  output$grid_plot_w1 <- renderPlotly({ plotly_w1() })
  
  
  ## heatmap W2 ----
  
  gg_plot_w2 <- reactive({
    req(input$use_window2_grid)
    res <- grid_results()
    req(res$w2)
    df <- copy(res$w2)
    
    frontier <- df[, .SD[which.max(richness_pct)], by = duration_h]
    
    dx_vals <- diff(sort(unique(df$start_off_h)))
    dy_vals <- diff(sort(unique(df$duration_h)))
    dx <- if (length(dx_vals) > 0) min(dx_vals) else 0.5
    dy <- if (length(dy_vals) > 0) min(dy_vals) else 0.5
    
    xmin <- min(df$start_off_h) - dx/2;  xmax <- max(df$start_off_h) + dx/2
    ymin <- min(df$duration_h)  - dy/2;  ymax <- max(df$duration_h)  + dy/2
    grid_v <- data.frame(x = unique(df$start_off_h) - dx/2)
    grid_h <- data.frame(y = unique(df$duration_h)  - dy/2)
    
    df[, tooltip := paste0(
      "<b>Window 2</b>",
      "<br>Start offset: ", round(start_off_h, 2), " h",
      "<br>Duration: ",     round(duration_h,  2), " h",
      "<br>Richness: ",     round(richness_pct, 1), " %",
      "<br>Species: ",      n_species_detected)]
    
    p <- ggplot(df, aes(x = start_off_h, y = duration_h)) +
      geom_contour(aes(z = richness_pct), binwidth = 5, color = "grey70", linewidth = 0.2) +
      geom_tile(aes(fill = richness_pct, text = tooltip), color = "grey89", linewidth = 0.5) +
      scale_fill_scico(palette = "bilbao", name = "Mean richness (%)", direction = -1) +
      labs(x = paste0("Start offset relative to ", input$grid_event2, " (hours)"),
           y = "Recording duration (hours)", title = "Window 2") +
      theme_minimal() +
      theme(panel.grid = element_blank(), axis.ticks = element_line(),
            axis.ticks.length = unit(3, "pt")) +
      scale_x_continuous(breaks = seq(floor(min(df$start_off_h)), ceiling(max(df$start_off_h)), by = 0.5)) +
      scale_y_continuous(breaks = seq(floor(min(df$duration_h)),  ceiling(max(df$duration_h)),  by = 0.5)) +
      geom_segment(data = grid_v, aes(x = x, xend = x, y = ymin, yend = ymax),
                   inherit.aes = FALSE, color = "grey89", linewidth = 0.1) +
      geom_segment(data = grid_h, aes(x = xmin, xend = xmax, y = y, yend = y),
                   inherit.aes = FALSE, color = "grey89", linewidth = 0.1)
    
    if (isTRUE(input$show_frontier)) {
      p <- p + geom_point(data = frontier,
                          aes(start_off_h, duration_h,
                              text = paste0("<b>Optimal frontier</b><br>",
                                            "Start offset: ", round(start_off_h, 2), " h<br>",
                                            "Duration: ", round(duration_h, 2), " h<br>",
                                            "Richness: ", round(richness_pct, 1), " %")),
                          inherit.aes = FALSE, color = "brown1", size = 3, shape = 8)
    }
    p
  })
  
  plotly_w2 <- reactive({ ggplotly(gg_plot_w2(), tooltip = "text") })
  output$grid_plot_w2 <- renderPlotly({ plotly_w2() })
  
  
  ## Scatter plot w1 x w2----
  
  gg_plot_comb <- reactive({
    res <- grid_results()
    req(res$w1, res$w2, res$comb)
    
    dt_raw_local <- copy(dt_raw())
    dt_raw_local[, date := as.Date(date)]
    n_recorder_days <- uniqueN(dt_raw_local[, .(recorder, date)])
    
    df1 <- copy(res$w1);   df1[, `:=`(type = "Window 1",  effort_total = duration_h * n_recorder_days * 60)]
    df2 <- copy(res$w2);   df2[, `:=`(type = "Window 2",  effort_total = duration_h * n_recorder_days * 60)]
    df3 <- copy(res$comb); df3[, `:=`(type = "Combined",  effort_total = effort * n_recorder_days)]
    
    df1[, tooltip := paste0("<b>Window 1</b><br>Start: ", round(start_off_h, 2), " h<br>",
                            "Duration: ", round(duration_h, 2), " h<br>",
                            "Effort: ", round(effort_total, 2), " h<br>",
                            "Richness: ", round(richness_pct, 1), " %<br>",
                            "Species: ", n_species_detected)]
    df2[, tooltip := paste0("<b>Window 2</b><br>Start: ", round(start_off_h, 2), " h<br>",
                            "Duration: ", round(duration_h, 2), " h<br>",
                            "Effort: ", round(effort_total, 2), " h<br>",
                            "Richness: ", round(richness_pct, 1), " %<br>",
                            "Species: ", n_species_detected)]
    df3[, tooltip := paste0("<b>Combined</b><br>Total effort: ", round(effort_total, 2), " h<br>",
                            "Richness: ", round(richness_pct, 1), " %<br>",
                            "Species: ", n_species_detected, "<br>",
                            "<b>W1</b>: start ", round(start_off_h, 2), " h, dur ", round(duration_h, 2), " h, richness ", round(richness_w1_pct, 1), "%<br>",
                            "<b>W2</b>: start ", round(start_off_h2, 2), " h, dur ", round(duration_h2, 2), " h, richness ", round(richness_w2_pct, 1), "%")]
    
    ggplot() +
      geom_point(data = df1, aes(x = effort_total, y = richness_pct, text = tooltip),
                 color = "steelblue", size = 2.5, alpha = 0.5) +
      geom_point(data = df2, aes(x = effort_total, y = richness_pct, text = tooltip),
                 color = "orange", size = 2.5, alpha = 0.5) +
      geom_point(data = df3, aes(x = effort_total, y = richness_pct, text = tooltip),
                 color = "darkgreen", size = 3.5, alpha = 0.7) +
      labs(x = "Total recording effort (minutes)", y = "Mean relative richness (%)") +
      theme_minimal()
  })
  
  plotly_comb <- reactive({ ggplotly(gg_plot_comb(), tooltip = "text") })
  output$grid_plot_combined <- renderPlotly({ plotly_comb() })
  
  ## plots organisation UI ----
  output$grid_plots_ui <- renderUI({
    if (isTRUE(input$use_window2_grid)) {
      div(style = "display: grid; grid-template-columns: 1fr; gap: 12px; width: 100%;",
          plotlyOutput("grid_plot_w1",       height = "450px", width = "100%"),
          plotlyOutput("grid_plot_w2",       height = "450px", width = "100%"),
          plotlyOutput("grid_plot_combined", height = "450px", width = "100%"))
    } else {
      plotlyOutput("grid_plot_w1", height = "700px", width = "100%")
    }
  })
  
  observeEvent(input$return_window, {
    assign("plotly_w1", plotly_w1(), envir = .GlobalEnv)
    if (isTRUE(input$use_window2_grid)) {
      assign("plotly_w2",   plotly_w2(),   envir = .GlobalEnv)
      assign("plotly_comb", plotly_comb(), envir = .GlobalEnv)
    }
  })
  
  observeEvent(input$return_window_data, {
    assign("time_window_data", grid_results(), envir = .GlobalEnv)
  })
  
  
  # TAB 4 MULTIOPTIMUM----
  
  multi_combo <- eventReactive(input$run_multi, {
    
    req(values$summary_richness, values$grid_results, values$schedules)
    message(Sys.time(), " multi_combo: start")
    
    summary_richness <- values$summary_richness
    schedules        <- values$schedules
    grid_results_raw <- values$grid_results
    two_windows      <- is.list(grid_results_raw) && !is.null(grid_results_raw$w2)
    
    window_sets <- if (two_windows) {
      list(W1 = grid_results_raw$w1, W2 = grid_results_raw$w2, Combined = grid_results_raw$comb)
    } else {
      list(W1 = if (is.list(grid_results_raw)) grid_results_raw$w1 else grid_results_raw)
    }
    
    dt      <- copy(dt_raw())
    req(dt)
    dt[, date := as.Date(date)]
    solar_dt <- setDT(solar_table())
    
    ## Pre-compute reference minute for W1 and W2 (done once, not per iteration) ----
    compute_ref_min <- function(event_input, solar) {
      if (event_input %in% c("dawn","sunrise","sunset","dusk")) {
        solar[, ref_min := as.integer(hour(get(event_input)) * 60 + minute(get(event_input)))]
      } else {
        fixed <- as.integer(sub(":.*", "", event_input)) * 60L +
          as.integer(sub(".*:", "", event_input))
        solar[, ref_min := fixed]
      }
      solar[, .(recorder, date, ref_min)]
    }
    
    ref_w1 <- compute_ref_min(input$grid_event,  copy(solar_dt))
    setnames(ref_w1, "ref_min", "ref_min_w1")
    
    if (two_windows && !is.null(input$grid_event2)) {
      ref_w2 <- compute_ref_min(input$grid_event2, copy(solar_dt))
      setnames(ref_w2, "ref_min", "ref_min_w2")
    } else {
      ref_w2 <- copy(ref_w1)
      setnames(ref_w2, "ref_min_w1", "ref_min_w2")
    }
    
    refs <- merge(ref_w1, ref_w2, by = c("recorder","date"), all = TRUE)
    
    ## Attach ref_min columns to each schedule (once per schedule) ----
    schedules_with_ref <- lapply(schedules, function(s) {
      merge(s, refs, by = c("recorder","date"), all.x = TRUE)
    })
    
    withProgress(message = "Running multi-optimum analysis...", value = 0, {
      
      combo_results <- list()
      k             <- 1L
      total_iter    <- sum(sapply(window_sets, function(gr) if (!is.null(gr)) nrow(gr) else 0L))
      done_iter     <- 0L
      
      for (set_name in names(window_sets)) {
        
        grid_r <- window_sets[[set_name]]
        if (is.null(grid_r)) next
        if (!is.data.table(grid_r)) setDT(grid_r)
        nW <- nrow(grid_r)
        if (nW == 0L) next
        if (set_name == "Combined") {
          top_n <- as.integer(input$top_n_combined)
          grid_r <- grid_r[order(duration_h, -richness_pct)][
            , .SD[seq_len(min(.N, top_n))], by = .(duration_h, duration_h2)
          ]
          nW <- nrow(grid_r)
        }
        for (w in seq_len(nW)) {
          
          done_iter <- done_iter + 1L
          if (done_iter %% 10L == 0L || done_iter == total_iter) {
            incProgress(10 / total_iter, detail = paste(done_iter, "/", total_iter))
          }
          
          win_row <- grid_r[w]
          
          ## Extract window offsets / durations by set type ----
          if (set_name == "W1") {
            w1_off_h <- win_row$start_off_h; w1_dur_h <- win_row$duration_h
            w2_off_h <- NA_real_;            w2_dur_h <- NA_real_
          } else if (set_name == "W2") {
            w1_off_h <- NA_real_;            w1_dur_h <- NA_real_
            w2_off_h <- win_row$start_off_h; w2_dur_h <- win_row$duration_h
          } else {
            w1_off_h <- win_row$start_off_h;  w1_dur_h <- win_row$duration_h
            w2_off_h <- win_row$start_off_h2; w2_dur_h <- win_row$duration_h2
          }
          
          for (p in seq_along(schedules_with_ref)) {
            
            sched <- copy(schedules_with_ref[[p]])
            
            ## Compute absolute window boundaries and filter blocks ----
            if (!is.na(w1_off_h)) sched[, `:=`(w1_start = ref_min_w1 + w1_off_h * 60,
                                               w1_end   = ref_min_w1 + (w1_off_h + w1_dur_h) * 60)]
            if (!is.na(w2_off_h)) sched[, `:=`(w2_start = ref_min_w2 + w2_off_h * 60,
                                               w2_end   = ref_min_w2 + (w2_off_h + w2_dur_h) * 60)]
            
            keep <- if (set_name == "W1") {
              sched[, block_in_window(start_min, end_min, w1_start, w1_end)]
            } else if (set_name == "W2") {
              sched[, block_in_window(start_min, end_min, w2_start, w2_end)]
            } else {
              sched[, block_in_window(start_min, end_min, w1_start, w1_end) |
                      block_in_window(start_min, end_min, w2_start, w2_end)]
            }
            
            sched_filt <- sched[keep]
            ## Keep only dates that are present in the detection data ----
            sched_filt <- sched_filt[as.Date(date) %in% unique(as.Date(dt$date))]
            
            ## Remove temporary boundary columns ----
            tmp_cols <- intersect(c("w1_start","w1_end","w2_start","w2_end","ref_min_w1","ref_min_w2"),
                                  names(sched_filt))
            if (length(tmp_cols)) sched_filt[, (tmp_cols) := NULL]
            
            if (nrow(sched_filt) == 0L) next
            
            res <- tryCatch(
              compute_richness(dt = dt, schedule = sched_filt,
                               per_spatial = "all", per_temporal = "day",
                               return_species = TRUE),
              error = function(e) NULL
            )
            
            richness_val       <- 0
            n_species_detected <- 0L
            
            if (!is.null(res) && !is.null(res$richness) && nrow(res$richness) > 0) {
              richness_val <- mean(res$richness$richness, na.rm = TRUE)
            }
            if (!is.null(res) && !is.null(res$species) && nrow(res$species) > 0) {
              sp_dt  <- as.data.table(res$species)
              sp_col <- if ("species" %in% names(sp_dt)) "species" else
                names(sp_dt)[which(sapply(sp_dt, function(x) is.character(x) || is.factor(x)))[1]]
             
            #   if (!is.null(sp_col) && "date" %in% names(sp_dt)) {
            #     n_species_detected <- round(
            #       mean(
            #         sp_dt[, .(n = uniqueN(get(sp_col))), by = date]$n
            #       ), 1)
            #   } else if (!is.null(sp_col)) {
            #     n_species_detected <- length(unique(sp_dt[[sp_col]]))
            #   }
            # }
              if (!is.null(sp_col)) n_species_detected <- length(unique(sp_dt[[sp_col]]))
            }
            
            combo_results[[k]] <- data.table(
              window_type        = set_name,
              start1_hour        = w1_off_h,
              dur1_h             = w1_dur_h,
              start2_hour        = w2_off_h,
              dur2_h             = w2_dur_h,
              duty_period        = attr(schedules[[p]], "period"),
              effort             = sum(sched_filt$end_min - sched_filt$start_min + 1L),
              richness           = richness_val,
              n_species_detected = n_species_detected
            )
            k <- k + 1L
          }
        }
      }
      
      if (length(combo_results) == 0L) {
        warning("No valid schedules produced any results.")
        return(NULL)
      }
      
      combo <- rbindlist(combo_results)
      combo[, richness_pct := 100 * richness / max(richness, na.rm = TRUE)]
      
      ## Compute Pareto front ----
      setorder(combo, effort)
      combo[, pareto := FALSE]
      best <- -Inf
      for (i in seq_len(nrow(combo))) {
        if (!is.na(combo$richness_pct[i]) && combo$richness_pct[i] > best) {
          combo$pareto[i] <- TRUE
          best <- combo$richness_pct[i]
        }
      }
      
      message(Sys.time(), " multi_combo: done")
      combo
    })
  })
  
  observeEvent(input$return_multi_data, {
    assign("multi_optimum_data", multi_combo(), envir = .GlobalEnv)
  })
  
  
  # TAB 4 PLOTS----
  
  gg_optimum <- reactive({
    dt <- multi_combo()
    req(dt)
    message(Sys.time(), " multi optimum: building plot")
    
    ## Recompute Pareto dominance for the plot ----
    dt[, dominated := FALSE]
    for (i in seq_len(nrow(dt))) {
      dt$dominated[i] <- any(
        dt$effort      <= dt$effort[i] &
          dt$richness_pct >= dt$richness_pct[i] &
          (dt$effort < dt$effort[i] | dt$richness_pct > dt$richness_pct[i])
      )
    }
    dt[, pareto := !dominated]
    pareto_front <- dt[pareto == TRUE][order(effort)]
    
    durations <- battery_specs$duration_h
    codes     <- battery_specs$code
    
    ###tooltips ----
    battery_tooltip <- function(eff_min) {
      paste0("<br><br><b>Full battery sets needed:</b><br>",
             paste0(codes, ": ", round((eff_min / 60) / durations, 2), collapse = "<br>"))
    }
    
    make_tooltip <- function(row, is_pareto = FALSE) {
      base <- paste0(
        if (is_pareto) "Local optimum<br>" else "",
        "Richness: ",        round(row$richness_pct, 1), "%<br>",
        "Species detected: ", row$n_species_detected, "<br>",
        "<b>Type:</b> ",      row$window_type, "<br>",
        "Effort: ",           row$effort, " min — ", round(row$effort / 60, 2), " h<br>"
      )
      w1_info <- if (!is.na(row$start1_hour))
        paste0("<b>Window 1</b><br>Start: ", round(row$start1_hour, 2), " h  Dur: ", round(row$dur1_h, 2), " h<br>")
      w2_info <- if (!is.na(row$start2_hour))
        paste0("<b>Window 2</b><br>Start: ", round(row$start2_hour, 2), " h  Dur: ", round(row$dur2_h, 2), " h<br>")
      paste0(base,
             if (!is.na(row$start1_hour)) w1_info else "",
             if (!is.na(row$start2_hour)) w2_info else "",
             "Period: ", row$duty_period,
             if (is_pareto) battery_tooltip(row$effort) else "")
    }
    
    dt[,           tooltip        := mapply(make_tooltip, split(dt,           seq_len(nrow(dt))),           FALSE)]
    pareto_front[, tooltip_pareto := mapply(make_tooltip, split(pareto_front, seq_len(nrow(pareto_front))), TRUE)]
    
    p <- ggplot(dt, aes(x = effort, y = richness_pct)) +
      geom_point(aes(color = window_type, shape = window_type, text = tooltip),
                 alpha = 0.7, size = 2) +
      geom_line(data = pareto_front, aes(group = 1), color = "red", linewidth = 0.5) +
      geom_point(data = pareto_front, aes(text = tooltip_pareto), color = "red", size = 2, alpha = 0.7) +
      scale_color_manual(values = c(W1 = "#1f78b4", W2 = "#6a3d9a", Combined = "#33a02c"),
                         name = "Window type") +
      scale_shape_manual(values = c(W1 = 16, W2 = 17, Combined = 18), name = "Window type") +
      labs(x = "Recording effort (min)", y = "Relative richness (%)",
           title = "Pareto frontier: optimal trade-off richness vs effort") +
      theme_minimal()
    
    p
  })
  
  plotly_optimum <- reactive({ ggplotly(gg_optimum(), tooltip = "text") })
  output$plot_multi_optimum <- renderPlotly({ plotly_optimum() })
  
  observeEvent(input$return_multi, {
    assign("plotly_optimum", plotly_optimum(), envir = .GlobalEnv)
  })
  
  
  # TAB 5  TARGET SPECIES ACTIVITY ----
  
  ## Populate species selector from current dataset ----
  observe({
    req(dt_raw())
    updateSelectizeInput(session, "target_species",
                         choices = sort(unique(dt_raw()$species)), server = TRUE)
  })
  
  target_result <- eventReactive(input$run_target, {
    
    req(dt_raw(), input$target_species)
    message(Sys.time(), " target species: start")
    
    dt   <- copy(dt_raw())
    recs <- unique(dt$recorder)
    solar_dt <- setDT(solar_table())
    
    sp <- dt[species %in% input$target_species]
    sp[, date := as.Date(date)]
    
    time_bin_min <- as.numeric(input$time_resolution)
    sp[, time_bin := (minute_of_day %/% time_bin_min) * time_bin_min]
    
    time_levels <- seq(0, 24 * 60 - time_bin_min, by = time_bin_min)
    
    ## Full grid of date × time_bin (ensures zeros are explicit) ----
    grid <- CJ(
      date     = seq(min(sp$date), max(sp$date), by = "day"),
      time_bin = time_levels
    )
    
    periods <- seq(input$target_period_min, input$target_period_max,
                   by = input$target_period_step)
    
    all_heatmaps <- list()
    
    withProgress(message = "Computing activity...", value = 0, {
      
      for (i in seq_along(periods)) {
        
        p <- periods[i]
        incProgress(1 / length(periods), detail = paste(i, "/", length(periods)))
        
        sched_rec_list <- lapply(recs, function(r) {
          sun_r  <- solar_dt[recorder == r]
          sched_r <- create_sun_duty_schedule(
            sun_table           = sun_r,
            period              = p,
            duty_duration       = input$target_duty_duration,
            nb_duty             = input$target_nb_duty,
            start_event         = "dawn",
            end_event           = "dusk",
            extend_before_hours = 1,
            extend_after_hours  = 1,
            record_24h          = TRUE
          )
          sched_r[, recorder := r]
          sched_r
        })
        
        sched <- rbindlist(sched_rec_list)
        
        ## Filter species detections to scheduled blocks ----
        sp_sched <- sched[sp, on = .(date, start_min <= minute_of_day, end_min >= minute_of_day),
                          nomatch = 0]
        
        counts <- sp_sched[, .(n = .N), by = .(date, time_bin)]
        counts <- merge(grid, counts, by = c("date","time_bin"), all.x = TRUE)
        counts[is.na(n), n := 0]
        
        counts[, `:=`(period_label = paste0("1/", p), period_num = p)]
        all_heatmaps[[as.character(p)]] <- counts
      }
    })
    
    heatmap_dt <- rbindlist(all_heatmaps)
    heatmap_dt[, period_label := factor(period_label,
                                        levels = paste0("1/", sort(unique(period_num))))]
    heatmap_dt[, date_chr    := as.character(date)]
    heatmap_dt[, hover_text  := paste0(
      "Time: ", sprintf("%02d:%02d", time_bin %/% 60, time_bin %% 60),
      "<br>Date: ", date_chr,
      "<br>Activity: ", n, " detections"
    )]
    
    message(Sys.time(), " target species: done")
    heatmap_dt
  })
  
  observeEvent(input$return_target_data, {
    assign("target_species_data", target_result(), envir = .GlobalEnv)
  })
  
  ## Animated plotly heatmap ----
  gg_target <- reactive({
    heatmap_dt <- target_result()
    req(heatmap_dt)
    
    global_max <- max(heatmap_dt$n, na.rm = TRUE)
    pal <- c("grey89", hcl.colors(1000, "Berlin"))
    
    plot_ly(
      data           = heatmap_dt,
      x              = ~date_chr,
      y              = ~time_bin / 60,
      z              = ~n,
      type           = "heatmap",
      frame          = ~period_label,
      colors         = pal,
      zmin           = 0,
      zmax           = global_max,
      text           = ~hover_text,
      hovertemplate  = "%{text}<extra></extra>",
      xgap           = 0.3,
      ygap           = 0.3,
      showscale      = TRUE,
      opacity        = 0.9
    ) %>%
      layout(
        title  = paste("Species detection:", paste(input$target_species, collapse = ", ")),
        yaxis  = list(title = "Hour of day", tickmode = "linear", tick0 = 0, dtick = 1,
                      range = c(0, 24), showgrid = FALSE),
        xaxis  = list(title = "Date", showgrid = FALSE, tickfont = list(size = 10), automargin = TRUE)
      ) %>%
      config(displaylogo = FALSE, modeBarButtonsToRemove = c("select2d","lasso2d"))
  })
  
  ## gg_target already returns a plotly object — render directly without ggplotly() ----
  output$target_heatmap <- renderPlotly({ gg_target() })
  
  observeEvent(input$return_target, {
    assign("target_species_plot", gg_target(), envir = .GlobalEnv)
  })
  
}


# Launch ----

shinyApp(ui, server)

