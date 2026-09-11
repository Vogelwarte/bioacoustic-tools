# FILTERING — BirdNET-ResChecker

Explore, filter, and visualize BirdNET detection data. Works with combined BirdNET selection tables.

With this app you can:
- Upload your BirdNET selection tables
- Apply custom confidence thresholds per species
- Visualize detections over time
- Export clean, filtered results (data and per-species summaries)

## Requirements

```r
install.packages(c(
  "shiny", "data.table", "readxl", "stringr", "plotly", "leaflet",
  "writexl", "RColorBrewer", "lubridate", "ggplot2", "scales",
  "bslib", "dplyr", "skimr", "DT"
))
```

## Usage

```r
shiny::runApp("FILTERING/1_App_1_Filtering.R")
```
or open `1_App_1_Filtering.R` in RStudio and click "Run App".

1. Upload one or more BirdNET selection tables (`.txt`).
2. If filenames contain coordinates (`Lat-XX_Long-XX_...`), check the GPS box to enable location-based features (map view).
3. Download the species XLS template, then fill in custom confidence thresholds per species and re-upload it — the graphs update automatically.
4. Optionally filter by date range or by recorder.
5. Use the sidebar buttons to download: raw counts, the filtered BirdNET file, and filtered per-species counts.

Tabs available: Overview, Visualizations, Filtering impact, Recorder Comparison, Recording Schedule, and a Map (when GPS mode is enabled).

---
Contributions: Amandine Serrurier, Jean-Nicolas Pradervand & Christophe Sahli — Swiss Ornithological Institute
