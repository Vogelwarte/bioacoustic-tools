# Bioacoustic tools

This repository hosts a suite of Shiny web applications for analyzing, filtering, visualizing, and optimizing bird detection results produced by the [BirdNET](https://github.com/kahst/BirdNET-Analyzer) algorithm. Together they cover the full workflow: preparing files, filtering and exploring detections, checking recording coverage, visualizing spectrograms, studying phenology, and optimizing recording schedules.

**Note:** these applications are independent third-party tools, developed and maintained separately from BirdNET. They are designed to work with BirdNET output files but are not affiliated with or endorsed by the BirdNET developers.

Works on Windows, macOS, and Linux.

## Contents

- [Apps overview](#apps-overview)
- [General requirements](#general-requirements)
- [RENAME](#rename--rename-files-with-coordinates-and-recorder-prefix)
- [FILTERING](#filtering--explore-and-filter-birdnet-detections)
- [RECTIME](#rectime--recording-time-coverage)
- [SPECTRO](#spectro--spectrograms--audio-playback)
- [PHENOLOGY](#phenology--temporal-activity-patterns)
- [SCORE](#score--sampling-coverage-optimizer-for-recording-effort)
- [Input files structure](#input-files-structure)
- [Common functions](#common-functions)
- [Contributors & acknowledgments](#contributors--acknowledgments)

## Apps overview

| App | What it does | Also available online |
|---|---|---|
| **RENAME** | Batch-renames audio/text files with a location + recorder + date/time prefix | — |
| **FILTERING** | Explores, filters (by confidence threshold), and visualizes BirdNET detections | [Filtering_Shiny_App](https://vogelwarte.shinyapps.io/Filtering_Shiny_App/) |
| **RECTIME** | Heatmap of recording coverage per day / site / recorder, to spot gaps | — |
| **SPECTRO** | Browses detections, generates spectrograms, plays matching audio | — |
| **PHENOLOGY** | Temporal activity graphs (phenology) with sunrise/sunset overlays | [PhenoApp_V1](https://vogelwarte.shinyapps.io/PhenoApp_V1/) |
| **SCORE** | Evaluates and optimizes recording schedules (duty cycle, time window) against species richness | [SCORE_online](https://vogelwarte.shinyapps.io/SCORE_online) |

Each app lives in its own folder (same name, all-caps) and can be run independently — you don't need to run them in order, except that **RENAME** is usually a good first step if your files don't yet follow the expected naming convention.

## General requirements

- R (version 4.1 or higher)
- RStudio (recommended)
- External dependencies for some apps: FFmpeg (MP3/FLAC support), Audacity (optional, for audio editing)

Each app lists its own required R packages below — most install automatically via `pacman` the first time the app runs.

---

## RENAME — Rename Files with Coordinates and Recorder Prefix

Batch-renames audio files by adding a location-based prefix (latitude/longitude, recorder, date, time). Processes every file in a selected directory at once. The prefix can be typed manually or generated automatically by clicking a location on the built-in map.

**⚠️ Required file name format:** `LatXX.XXXX_LongYY.YYYY_Recorder_Date_Time.wav`

**Requirements**
```r
install.packages(c("shiny", "shinyFiles", "leaflet", "bslib"))
```

**Usage**
```r
shiny::runApp("RENAME/0_App_0_Rename.R")
```
1. Click **"Choose Folder"** to select the directory containing the files to rename.
2. Get a prefix either by clicking on the map (auto lat/long) or by typing it manually.
3. Click **"Preview Changes"** to check the new file names before applying anything.
4. Click **"Apply Renaming"**.

**⚠️ Warning:** renaming is irreversible — original file names cannot be restored afterwards.

---

## FILTERING — Explore and filter BirdNET detections

Explores, filters, and visualizes BirdNET detection data. Works with combined BirdNET selection tables.

With this app you can: upload your BirdNET selection tables, apply custom confidence thresholds per species, visualize detections over time, and export clean, filtered results.

Also available online, no R install needed: **[Filtering_Shiny_App](https://vogelwarte.shinyapps.io/Filtering_Shiny_App/)** — note the online version is considerably slower than running it locally.

**Requirements**
```r
install.packages(c(
  "shiny", "data.table", "readxl", "stringr", "plotly", "leaflet",
  "writexl", "RColorBrewer", "lubridate", "ggplot2", "scales",
  "bslib", "dplyr", "skimr", "DT"
))
```

**Usage**
```r
shiny::runApp("FILTERING/1_App_1_Filtering.R")
```
1. Upload one or more BirdNET selection tables (`.txt`).
2. If filenames contain coordinates (`Lat-XX_Long-XX_...`), check the GPS box to enable location-based features (map view).
3. Download the species XLS template, fill in custom confidence thresholds per species, and re-upload it — graphs update automatically.
4. Optionally filter by date range or by recorder.
5. Download results (raw counts, filtered BirdNET file, filtered per-species counts) from the sidebar.

Tabs available: Overview, Visualizations, Filtering impact, Recorder Comparison, Recording Schedule, and a Map (when GPS mode is enabled).

---

## RECTIME — Recording time coverage

Scans a local folder of BirdNET selection tables (or audio files) and shows a heatmap of recording coverage by day, per site/recorder — useful for spotting missing data.

**Requirements**
```r
install.packages(c("shiny", "ggplot2", "bslib", "rstudioapi"))
```

**Usage**

Open `4_app_4_recording_time.R` in RStudio and click **Run App**.

1. Paste a folder path (or use **Browse**), then click **Scan folder**.
2. R lists the files, reads the date/time from each filename, and infers site/recorder from the folder structure.
3. Use the filters (file type, site, date range) to narrow things down.
4. The heatmap shows one row per recorder — green = day with a recording, grey = a gap.
5. Download the plot (PNG) or the data (CSV) from the sidebar.

**Note:** everything runs locally — R needs direct access to the folder on disk. This won't work if the app is hosted on a server for other people to use from their own computers.

---

## SPECTRO — Spectrograms & audio playback

Analyzes BirdNET detection results: browse detections, generate spectrograms, and listen to the corresponding audio.

**Requirements**
```r
install.packages(
  "shiny", "grid", "bslib", "stringr", "data.table", "shinyFiles", "tidyr",
  "stringi", "fs", "parallel", "ggplot2", "mgcv", "gratia", "hms", "suncalc",
  "lubridate", "DT", "seewave", "tuneR", "base64enc", "av", "audio",
  "future.apply", "signal"
)
```
Needs the `Common_functions` folder alongside the app to load its helper functions.

**Usage**
```r
shiny::runApp("SPECTRO/2_App_2_Spectro.R")
```
1. Select the BirdNET results folder and the corresponding audio folder.
2. Set the device and deployment timezones, then click **"Start App"**.
3. In the "Spectrogram & Export" tab, adjust the High-Pass filter and frequency range to explore a detection's spectrogram and play the matching audio clip.

---

## PHENOLOGY — Temporal activity patterns

Analyzes BirdNET detection results and produces phenology graphs (detection activity across the season / day), with sunrise/sunset overlays.

Also available online, no R install needed: **[PhenoApp_V1](https://vogelwarte.shinyapps.io/PhenoApp_V1/)** — note the online version is considerably slower than running it locally.

**Requirements**
```r
install.packages(c(
  "shiny", "grid", "bslib", "stringr", "data.table", "shinyFiles", "tidyr",
  "janitor", "stringi", "fs", "parallel", "ggplot2", "mgcv", "gratia", "hms",
  "suncalc", "scales", "lubridate", "DT", "future.apply"
))
```
Needs the `Common_functions` folder alongside the app to load its helper functions.

**Usage**
```r
shiny::runApp("PHENOLOGY/3_App_3_Phenology.R")
```
1. Select the BirdNET results folder (check "Results in compiled format?" if relevant), then click **"Start App"**.
2. Set the device and deployment timezones and click **"Update Timezones"** if they differ.
3. In the Phenology panel: pick a species, set the confidence threshold, latitude/longitude, and aggregation interval (minutes).
4. Click **"Generate Plot"**.

---

## SCORE — Sampling Coverage Optimizer for Recording Effort

Evaluates and optimizes passive acoustic sampling designs from an existing BirdNET dataset. It simulates how different recording schedules (duty cycle, time window, or both) would perform on the same data, to help decide on recording protocols while retaining a high proportion of species richness.

Intended for a pilot BirdNET dataset of a few days or more, for one or several recorders. Three example datasets are included in the `Example datasets` folder.

Also available online, no R install needed: **[SCORE_online](https://vogelwarte.shinyapps.io/SCORE_online)** — the online version is not recommended for real analyses, as computation time is strongly affected and input file size may be limited.

**Requirements**
```r
install.packages(c(
  "shiny", "shinythemes", "data.table", "lubridate",
  "ggplot2", "suncalc", "dplyr", "DT", "scico",
  "plotly", "stringr", "readxl"
))
```

**Usage**
```r
shiny::runApp("path/to/SCORE.R")
```
General workflow: import and prepare the BirdNET dataset → define coordinates/timezone → apply filters → explore duty-cycle designs → explore recording-window designs → compare complete designs (multi-optimum) → optionally visualize target-species activity under different schedules.

SCORE is the most extensive app of the suite, with a detailed tab-by-tab guide (including screenshots) already available directly in **`SCORE/README.md`**, and a full **User Guide** tab inside the app itself — see those for the complete documentation, this section is only a summary.

---

## Input files structure

Expected input files:
- Detection results: `*.BirdNET.selection.table.txt` or compiled selection tables: `*BirdNET_SelectionTable.txt`
- Audio files: `*.wav`, `*.WAV`, `*.mp3`, `*.flac` (in matching folder structure)

## Common functions

Shared helper functions used by several apps are located in the `Common_functions` folder. They are sourced automatically when running an app, as long as `Common_functions` is in the same location as the app's `.R` file.

## Online Apps

Some of the Shiny Apps are available online. You can access these different Apps via le links below: 
- Filtering App: https://vogelwarte.shinyapps.io/Filtering_Shiny_App/
- Phenology App: https://vogelwarte.shinyapps.io/PhenoApp_V1/
- SCORE App:     

## Contributors & acknowledgments

Contributions: Amandine Serrurier, Jean-Nicolas Pradervand & Christophe Sahli
MIT License © 2026 vogelwarte.ch

This application is designed to work with output from [BirdNET](https://github.com/kahst/BirdNET-Analyzer), a powerful AI-based bird sound identification system developed by the K. Lisa Yang Center for Conservation Bioacoustics at the Cornell Lab of Ornithology.
