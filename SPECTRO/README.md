# SPECTRO TOOL APP

Analyzes BirdNET detection results: browse detections, generate spectrograms, and listen to the corresponding audio.

## Requirements

```r
install.packages(
  "shiny", "grid", "bslib", "stringr", "data.table", "shinyFiles", "tidyr",
  "stringi", "fs", "parallel", "ggplot2", "mgcv", "gratia", "hms", "suncalc",
  "lubridate", "DT", "seewave", "tuneR", "base64enc", "av", "audio",
  "future.apply", "signal"
)
```


Needs the `Common_functions` folder alongside the app to load its helper functions.

## Usage

```r
shiny::runApp("SPECTRO/2_App_2_Spectro.R")
```
or open `2_App_2_Spectro.R` in RStudio and click "Run App".

1. Select the BirdNET results folder and the corresponding audio folder.
2. Set the device and deployment timezones, then click **"Start App"**.
3. In the "Spectrogram & Export" tab, adjust the High-Pass filter and the frequency range to explore a detection's spectrogram and play the matching audio clip.

---

