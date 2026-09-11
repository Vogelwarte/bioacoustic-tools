# Recording Time

A Shiny app that scans a local folder of BirdNET selection tables (or audio
files) and shows a heatmap of recording coverage by day, per site/recorder.

## Packages needed

Install once:

```r
install.packages(c("shiny", "ggplot2", "bslib"))
```

Optional, only for the "Browse" folder-picker button (RStudio only):

```r
install.packages("rstudioapi")
```

## How to run

Open `app.R` in RStudio and click **Run App**.

## How it works

1. Paste a folder path (or use **Browse**), then click **Scan folder**.
2. R lists the files in that folder, reads the date/time in each filename,
   and figures out the site and recorder from the folder structure.
3. Use the filters (file type, site, date range) to narrow things down.
4. The heatmap shows one row per recorder: green = day with a recording,
   grey = a gap.
5. Download the plot (PNG) or the data (CSV) from the sidebar.

**Note:** everything runs locally — R needs direct access to the folder on
disk. This won't work if the app is hosted on a server for other people to
use from their own computers.
