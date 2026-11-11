save_run <- function(run_id, threshold = 90,
                     dir = "/Users/amandine.serrurier/Documents/SCORE_runs",
                     clean = TRUE) {
  library(data.table)
  d <- file.path(dir, run_id)
  dir.create(file.path(d, "plots"),  recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(d, "tables"), recursive = TRUE, showWarnings = FALSE)
  
  ## 1. Paramètres lisibles ----
  p <- multi_params
  keep <- c("timezone", "min_conf", "record_24h", "period_min", "period_max", "period_step",
            "duty_duration", "grid_event", "grid_start_min", "grid_start_max", "grid_start_step",
            "grid_dur_min", "grid_dur_max", "grid_dur_step", "use_window2_grid", "grid_event2",
            "grid_start_min2", "grid_start_max2", "grid_start_step2", "grid_dur_min2",
            "grid_dur_max2", "grid_dur_step2", "top_n_combined")
  params <- data.table(parameter = keep,
                       value = sapply(keep, function(k) paste(p[[k]], collapse = ", ")))
  params <- rbind(params, data.table(parameter = c("first_date", "last_date", "n_days"),
                                     value = c(as.character(multi_dates), diff(multi_dates) + 1)))
  fwrite(params, file.path(d, "parameters.csv"))
  
  ## 2. Meilleur schedule ----
  res  <- as.data.table(multi_optimum_data)
  n_days <- as.integer(diff(multi_dates)) + 1
  n_rec  <- if (exists("multi_recorders")) length(multi_recorders) else 1
  res[, effort_pct_24h       := 100 * effort / (n_days * 1440 * n_rec)]
  res[, effort_reduction_pct := 100 - effort_pct_24h]
  best <- res[richness_pct >= threshold][order(effort)][1]
  fwrite(res,  file.path(d, "tables", "multi_optimum_all_candidates.csv"))
  fwrite(best, file.path(d, "best_schedule.csv"))
  
  ## 3. Tous les graphes disponibles ----
  plots <- c("ggplot_rarefaction", "gg_heatmap_missed",                     # onglet Duty cycle
             "ggplot_w1", "ggplot_w2", "ggplot_comb",                       # onglet Time window
             "ggplot_optimum",                                              # onglet Multi-optimum
             "plotly_rarefaction", "plotly_heatmap_missed",
             "plotly_w1", "plotly_w2", "plotly_comb", "plotly_optimum",
             "target_species_plot")                                         # onglet Target (plotly)
  saved <- c()
  for (nm in plots) {
    if (!exists(nm, envir = .GlobalEnv)) next
    obj <- get(nm, envir = .GlobalEnv)
    saveRDS(obj, file.path(d, "plots", paste0(nm, ".rds")))      # objet R, ré-ouvrable
    if (inherits(obj, "ggplot")) {
      ggplot2::ggsave(file.path(d, "plots", paste0(nm, ".png")), obj,
                      width = 8, height = 6, dpi = 300)
    } else if (inherits(obj, "plotly")) {
      htmlwidgets::saveWidget(obj, file.path(d, "plots", paste0(nm, ".html")),
                              selfcontained = TRUE)              # interactif, s'ouvre dans le navigateur
    }
    saved <- c(saved, nm)
  }
  
  ## 4. Autres tables de résultats ----
  for (nm in c("duty_cycle_data", "time_window_data", "target_species_data")) {
    if (exists(nm, envir = .GlobalEnv))
      saveRDS(get(nm, envir = .GlobalEnv), file.path(d, "tables", paste0(nm, ".rds")))
  }
  
  ## 5. Ligne de résumé dans l'index ----
  idx <- data.table(run_id, first_date = multi_dates[1], last_date = multi_dates[2],
                    n_days = as.integer(diff(multi_dates)) + 1, threshold, best,
                    plots_saved = paste(saved, collapse = "; "))
  f <- file.path(dir, "index.csv")
  fwrite(idx, f, append = file.exists(f))
  
  ## 6. Nettoyage : évite de ré-enregistrer les graphes de ce run avec le suivant ----
  if (clean) {
    rm(list = intersect(c(plots, "multi_params", "multi_dates", "multi_optimum_data",
                          "duty_cycle_data", "time_window_data", "target_species_data"),
                        ls(.GlobalEnv)), envir = .GlobalEnv)
  }
  message("Enregistré dans : ", d, "\nGraphes : ", paste(saved, collapse = ", "))
  invisible(best)
}
save_run("0106-0706_7day")

