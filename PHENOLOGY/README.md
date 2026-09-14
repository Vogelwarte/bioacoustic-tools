# PHENOLOGY APP TOOL

Analyzes BirdNET detection results and produces phenology graphs (detection activity across the season / day).

## Requirements

```r
install.packages(
  "shiny", "grid", "bslib", "stringr", "data.table", "shinyFiles", "tidyr",
  "janitor", "stringi", "fs", "parallel", "ggplot2", "mgcv", "gratia", "hms",
  "suncalc", "scales", "lubridate", "DT", "future.apply"
)
```


Needs the `Common_functions` folder alongside the app to load its helper functions.

## Usage

```r
shiny::runApp("PHENOLOGY/3_App_3_Phenology.R")
```
or open `3_App_3_Phenology.R` in RStudio and click "Run App".

1. Select the BirdNET results folder (check "Results in compiled format?" if relevant), then click **"Start App"**.
2. Set the device and deployment timezones and click **"Update Timezones"** if they differ.
3. In the Phenology panel: pick a species, set the confidence threshold, latitude/longitude, and aggregation interval (minutes).
4. Click **"Generate Plot"** to produce the phenology graph.

---
