# Recording Time


library(shiny)
library(ggplot2)

# ---- Detection site / recorder from path -------------------
# Exemple : .../traning_data_complete_bn_2024_GG/Gondo/GON_A/1/Data/SMA14439_...


is_generic_folder <- function(x) {
  y <- tolower(trimws(x))
  generic <- c(
    "data", "audio", "recording", "recordings", "result", "results",
    "output", "outputs", "birdnet", "files", "volumes", "transfer",
    "analysis_results", "analisys_results"
  )
  y %in% generic ||
    grepl("^(sd|card)?[ _-]*[0-9]+$", y) ||
    grepl("^20[0-9]{2}$", y) ||
    grepl("^(visit|session|deployment|round|run)[ _-]*[0-9]+$", y)
}

infer_site_recorder <- function(dir_path) {
  parts <- strsplit(dir_path, "/", fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]

  recorder_pattern <- "^[A-Za-z0-9-]+_[A-Za-z0-9]{1,3}$"
  matches <- which(grepl(recorder_pattern, parts))
  rec_pos <- if (length(matches)) max(matches) else 0L

  site <- "Unknown"
  recorder <- "Unknown"

  if (rec_pos > 0) {
    recorder <- parts[rec_pos]
    site <- sub("_[^_]+$", "", recorder)

    if (rec_pos > 1) {
      for (i in (rec_pos - 1):1) {
        if (!is_generic_folder(parts[i])) {
          site <- parts[i]
          break
        }
      }
    }
  } else {
    useful <- Filter(Negate(is_generic_folder), parts)
    n <- length(useful)
    if (n >= 2) {
      site <- useful[n - 1]
      recorder <- useful[n]
    } else if (n == 1) {
      recorder <- useful[1]
      site <- sub("_[^_]+$", "", recorder)
    }
  }

  c(site = site, recorder = recorder)
}

# ---- Scan du dossier ---------------------------------------------------------

matches_file_type <- function(name, type) {
  if (type == "audio")   return(grepl("\\.(wav|flac|mp3|m4a|aif|aiff)$", name, ignore.case = TRUE))
  if (type == "birdnet") return(grepl("selection\\.table\\.txt$", name, ignore.case = TRUE))
  TRUE
}

# Cette fonction fait tout le travail : lister les fichiers, ne garder que
# ceux du bon type avec un horodatage dans leur nom, en déduire site/recorder/
# date, puis regrouper par jour. incProgress() alimente la barre de
# progression native de Shiny (appelée depuis withProgress() dans le server).
scan_folder <- function(root, file_type) {
  rel <- list.files(root, recursive = TRUE, full.names = FALSE)
  incProgress(1 / 3, detail = paste(length(rel), "files listed"))

  name <- basename(rel)
  keep <- !grepl("\\.zip$", name, ignore.case = TRUE) & matches_file_type(name, file_type)
  total <- length(keep)
  rel <- rel[keep]
  name <- name[keep]

  # Horodatage AAAAMMJJ[_-]HHMMSS (ou HHMM) dans le nom du fichier
  stamp_pattern <- "20[0-9]{6}[_-]?([0-9]{6}|[0-9]{4})"
  has_stamp <- grepl(stamp_pattern, name)
  rel <- rel[has_stamp]
  name <- name[has_stamp]

  incProgress(1 / 3, detail = paste(length(rel), "matching files"))

  if (!length(rel)) {
    return(list(total = total, matched = 0L, daily = data.frame(), preview = data.frame()))
  }

  stamp <- gsub("[^0-9]", "", regmatches(name, regexpr(stamp_pattern, name)))
  stamp <- ifelse(nchar(stamp) == 12, paste0(stamp, "00"), stamp)
  date_str <- substr(stamp, 1, 8)

  # Site / recorder calculés une seule fois par dossier (pas par fichier),
  # puis reportés sur chaque fichier de ce dossier.
  dirs <- dirname(rel)
  uniq_dirs <- unique(dirs)
  meta <- vapply(uniq_dirs, infer_site_recorder, c(site = "", recorder = ""))
  idx <- match(dirs, uniq_dirs)
  site <- gsub("[\t\r\n]", " ", meta["site", idx])
  recorder <- gsub("[\t\r\n]", " ", meta["recorder", idx])

  incProgress(1 / 3, detail = "Grouping by day...")

  daily <- aggregate(
    list(recordings = rep(1L, length(rel))),
    by = list(site = site, recorder = recorder, date = date_str),
    FUN = sum
  )

  n_preview <- min(12L, length(rel))
  preview <- data.frame(
    site = site[seq_len(n_preview)],
    recorder = recorder[seq_len(n_preview)],
    timestamp = stamp[seq_len(n_preview)],
    file = name[seq_len(n_preview)],
    stringsAsFactors = FALSE
  )

  list(total = total, matched = length(rel), daily = daily, preview = preview)
}

# ---- Helpers UI ---------------------------------------------------------------

stat_box <- function(value, label, width) {
  column(
    width,
    wellPanel(
      tags$h3(value, style = "margin-top:0; margin-bottom:4px;"),
      tags$span(label)
    )
  )
}

# ---- UI -----------------------------------------------------------------------

ui <- fluidPage(
  theme = bslib::bs_theme(version = 3, bootswatch = "darkly"),
  div(
    style = "display: flex; justify-content: space-between; align-items: center; padding: 10px 20px; background-color: #0f0f0f; border-bottom: 1px solid #333; margin-bottom: 20px;",
    div(style = "font-size: 24px; font-weight: bold; color: #fff;", "Bioacoustic tools - Recording Time"),
    div(tags$img(src = "logo.png", height = "70px", style = "margin-right: 10px; border-radius: 12px; padding: 6px 10px; background-color: white;"))
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      tags$label("Folder to scan", `for` = "folder_path"),
      fluidRow(
        column(9, textInput("folder_path", NULL, placeholder = "/Users/.../Data")),
        column(3, actionButton("browse", NULL, icon = icon("folder-open"), title = "Browse (RStudio only)"))
      ),
      actionButton("scan", "Scan folder", icon = icon("search"), class = "btn-primary", width = "100%"),
      tags$p(textOutput("scan_status"), class = "text-muted", style = "margin-top:8px; min-height:18px;"),

      selectInput(
        "file_type",
        "File type",
        choices = c(
          "BirdNET selection tables" = "birdnet",
          "Audio files" = "audio",
          "Any timestamped file" = "any"
        ),
        selected = "birdnet"
      ),

      selectInput("site_filter", "Site", choices = "All sites"),
      uiOutput("date_range_ui"),

      radioButtons(
        "metric",
        "Heatmap value",
        choices = c(
          "Number of recordings per day" = "files",
          "Estimated recording hours per day" = "hours"
        ),
        selected = "files"
      ),

      conditionalPanel(
        "input.metric == 'hours'",
        numericInput("minutes_per_file", "Minutes per file", value = 60, min = 0.01, step = 1)
      ),

      hr(),
      downloadButton("download_plot", "Download PNG"),
      downloadButton("download_csv", "Download CSV")
    ),

    mainPanel(
      width = 9,
      uiOutput("summary_boxes"),
      plotOutput("coverage_plot", height = "650px"),
      tags$details(
        tags$summary("File interpretation preview"),
        br(),
        tableOutput("preview")
      )
    )
  )
)

# ---- Server ---------------------------------------------------------------

server <- function(input, output, session) {

  # ---- Choix du dossier -------------------------------------------------

  observeEvent(input$browse, {
    if (requireNamespace("rstudioapi", quietly = TRUE) && isTRUE(rstudioapi::isAvailable())) {
      chosen <- tryCatch(rstudioapi::selectDirectory(caption = "Choose a folder"), error = function(e) NULL)
      if (!is.null(chosen)) updateTextInput(session, "folder_path", value = chosen)
    } else {
      showNotification("Folder picker unavailable here — paste the path directly in the text field.", type = "warning")
    }
  })

  # ---- Scan (uniquement au clic sur "Scan folder") -----------------------
  # eventReactive : contrairement à reactive(), ne se relance QUE sur cet
  # événement précis (le clic du bouton), jamais tout seul.

  folder_data <- eventReactive(input$scan, {
    path <- input$folder_path
    validate(need(nzchar(path) && dir.exists(path), "Folder not found — check the path above."))

    raw <- withProgress(message = "Scanning folder...", value = 0, {
      scan_folder(path, input$file_type)
    })

    daily <- raw$daily
    if (nrow(daily)) {
      daily$date <- as.Date(daily$date, "%Y%m%d")

      # Si un nom de site a été récupéré par erreur comme "recorder" sous un
      # autre site, on retire cette association croisée (les cas légitimes où
      # site == recorder restent valides).
      site_names <- unique(daily$site)
      wrong_site_as_recorder <- daily$recorder %in% site_names & daily$recorder != daily$site
      daily <- daily[!wrong_site_as_recorder, , drop = FALSE]
    }

    preview <- raw$preview
    if (nrow(preview)) {
      preview$timestamp <- as.POSIXct(preview$timestamp, "%Y%m%d%H%M%S", tz = "UTC")
      if (nrow(daily)) {
        valid_pairs <- unique(paste(daily$site, daily$recorder, sep = "\r"))
        preview_pairs <- paste(preview$site, preview$recorder, sep = "\r")
        preview <- preview[preview_pairs %in% valid_pairs, , drop = FALSE]
      }
    }

    list(
      daily = daily,
      preview = preview,
      total = raw$total,
      matched = if (nrow(daily)) sum(daily$recordings) else 0L
    )
  }, ignoreInit = TRUE)

  output$scan_status <- renderText({
    info <- folder_data()
    paste0(format(info$matched, big.mark = ","), " / ", format(info$total, big.mark = ","), " matching files")
  })

  observeEvent(folder_data(), {
    dat <- folder_data()$daily

    if (!nrow(dat)) {
      updateSelectInput(session, "site_filter", choices = "All sites")
      return()
    }

    updateSelectInput(
      session,
      "site_filter",
      choices = c("All sites", sort(unique(dat$site))),
      selected = "All sites"
    )
  }, ignoreInit = TRUE)

  # Le sélecteur de dates est régénéré à chaque nouveau scan, avec min/max
  # calés exactement sur les dates trouvées : c'est ce qui garantit qu'il
  # reste toujours aligné avec ce que montre le graphique. (On évite
  # updateDateRangeInput(), dont la mise à jour de min/max n'est pas fiable
  # à 100% une fois le widget déjà affiché.)
  output$date_range_ui <- renderUI({
    dat <- folder_data()$daily
    req(nrow(dat) > 0)

    full <- range(dat$date, na.rm = TRUE)

    dateRangeInput(
      "dates",
      "Date range",
      start = full[1],
      end = full[2],
      min = full[1],
      max = full[2],
      separator = " to "
    )
  })

  selected_data <- reactive({
    dat <- folder_data()$daily
    req(nrow(dat) > 0)

    if (!is.null(input$site_filter) && input$site_filter != "All sites") {
      dat <- dat[dat$site == input$site_filter, , drop = FALSE]
    }

    dat
  })

  selected_dates <- reactive({
    dat <- selected_data()
    full <- range(dat$date, na.rm = TRUE)
    chosen <- as.Date(input$dates)

    if (length(chosen) != 2 || anyNA(chosen)) {
      return(full)
    }

    start <- max(chosen[1], full[1])
    end <- min(chosen[2], full[2])

    # Peut arriver une fraction de seconde pendant le chargement d'un nouveau dossier.
    if (start > end) {
      return(full)
    }

    c(start, end)
  })

  heat_data <- reactive({
    dat <- selected_data()
    r <- selected_dates()

    dat <- dat[dat$date >= r[1] & dat$date <= r[2], , drop = FALSE]
    req(nrow(dat) > 0)

    recorders <- unique(dat[c("site", "recorder")])
    dates <- data.frame(date = seq(r[1], r[2], by = "day"))

    grid <- merge(recorders, dates, by = NULL)
    grid <- merge(grid, dat, by = c("site", "recorder", "date"), all.x = TRUE, sort = FALSE)
    grid$recordings[is.na(grid$recordings)] <- 0L

    grid$value <- if (input$metric == "hours") {
      grid$recordings * input$minutes_per_file / 60
    } else {
      grid$recordings
    }

    grid
  })

  coverage_plot <- reactive({
    dat <- heat_data()
    req(nrow(dat) > 0)

    r <- selected_dates()
    span <- as.integer(r[2] - r[1])
    date_breaks <- if (span <= 45) {
      "1 week"
    } else if (span <= 200) {
      "2 weeks"
    } else if (span <= 450) {
      "1 month"
    } else {
      "2 months"
    }
    date_labels <- if (format(r[1], "%Y") == format(r[2], "%Y")) "%d %b" else "%d %b\n%Y"

    positive <- dat[dat$recordings > 0, , drop = FALSE]
    zero <- dat[dat$recordings == 0, , drop = FALSE]

    max_value <- if (nrow(positive)) max(positive$value, na.rm = TRUE) else 1
    if (!is.finite(max_value) || max_value <= 0) max_value <- 1

    legend_title <- if (input$metric == "hours") "Hours/day" else "Recordings/day"
    status_labels <- c("Green = day with recording", "Grey = day without recording")

    # Points invisibles utilisés uniquement pour créer la légende à 2 items.
    # "site" est indispensable ici : sans elle, facet_grid ne sait pas à quel
    # panneau (site) rattacher ces points invisibles, et les dessine dans TOUS
    # les panneaux — créant une ligne "fantôme" du nom du recorder pour
    # chaque site. On réutilise une paire (site, recorder) déjà réelle dans
    # les données, donc aucune nouvelle catégorie n'est ajoutée nulle part.
    legend_key <- data.frame(
      site = rep(dat$site[1], 2),
      date = rep(r[1], 2),
      recorder = rep(dat$recorder[1], 2),
      status = factor(status_labels, levels = status_labels)
    )

    ggplot() +
      geom_tile(
        data = zero,
        aes(date, recorder),
        fill = "#D9D9D6",
        colour = "#FFFFFF",
        linewidth = 0.18,
        width = 1,
        height = 0.82
      ) +
      geom_tile(
        data = positive,
        aes(date, recorder, fill = value),
        colour = "#FFFFFF",
        linewidth = 0.18,
        width = 1,
        height = 0.82
      ) +
      geom_point(
        data = legend_key,
        aes(date, recorder, colour = status),
        alpha = 0,
        size = 0,
        inherit.aes = FALSE,
        show.legend = TRUE
      ) +
      scale_fill_gradient(
        low = "#C6D7C8",
        high = "#52745B",
        limits = c(0, max_value),
        name = legend_title
      ) +
      scale_colour_manual(
        values = c(
          "Green = day with recording" = "#6F9278",
          "Grey = day without recording" = "#D9D9D6"
        ),
        breaks = status_labels,
        name = NULL
      ) +
      scale_x_date(
        date_breaks = date_breaks,
        date_labels = date_labels,
        expand = c(0, 0)
      ) +
      facet_grid(site ~ ., scales = "free_y", space = "free_y", switch = "y") +
      labs(
        title = "Recording coverage",
        subtitle = paste0(
          format(r[1], "%d %b %Y"), " to ", format(r[2], "%d %b %Y"),
          "  •  ", format(sum(dat$recordings), big.mark = ","), " recordings"
        ),
        x = "Date",
        y = "Recorder"
      ) +
      guides(
        fill = guide_colourbar(order = 1),
        colour = guide_legend(
          order = 2,
          override.aes = list(shape = 15, size = 5, alpha = 1)
        )
      ) +
      theme_minimal(base_size = 13) +
      theme(
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        panel.spacing.y = grid::unit(0.7, "lines"),
        plot.title = element_text(face = "bold")
      )
  })

  output$summary_boxes <- renderUI({
    dat <- heat_data()
    r <- selected_dates()

    n_files <- sum(dat$recordings)
    minutes <- if (is.null(input$minutes_per_file)) 60 else input$minutes_per_file

    fluidRow(
      stat_box(format(n_files, big.mark = ","), "Files", 2),
      stat_box(length(unique(dat$site)), "Sites", 2),
      stat_box(length(unique(dat$recorder)), "Recorders", 2),
      stat_box(paste0(as.integer(r[2] - r[1]) + 1, " days"), "Period", 3),
      stat_box(format(round(n_files * minutes / 60, 1), big.mark = ","), "Estimated hours", 3)
    )
  })

  output$coverage_plot <- renderPlot({
    coverage_plot()
  }, res = 120)

  output$preview <- renderTable({
    dat <- folder_data()$preview
    req(nrow(dat) > 0)

    data.frame(
      Site = dat$site,
      Recorder = dat$recorder,
      Timestamp = format(dat$timestamp, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
      File = dat$file,
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE, rownames = FALSE)

  output$download_plot <- downloadHandler(
    filename = function() paste0("recording_time_", Sys.Date(), ".png"),
    content = function(file) {
      dat <- heat_data()
      height <- max(6, min(30, 2.5 + 0.28 * length(unique(dat$recorder))))
      ggsave(file, plot = coverage_plot(), width = 14, height = height, dpi = 220, bg = "white")
    }
  )

  output$download_csv <- downloadHandler(
    filename = function() paste0("recording_time_", Sys.Date(), ".csv"),
    content = function(file) {
      dat <- heat_data()
      minutes <- if (is.null(input$minutes_per_file)) 60 else input$minutes_per_file

      write.csv(
        data.frame(
          site = dat$site,
          recorder = dat$recorder,
          date = dat$date,
          recordings = dat$recordings,
          estimated_hours = round(dat$recordings * minutes / 60, 4),
          has_recording = dat$recordings > 0
        ),
        file,
        row.names = FALSE
      )
    }
  )
}

shinyApp(ui, server)
