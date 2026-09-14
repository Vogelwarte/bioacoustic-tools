# RENAME APP TOOL

Batch-renames audio files by adding a location-based prefix (latitude/longitude, recorder, date, time). Processes every file in a selected directory at once.

The prefix can be typed manually or generated automatically by clicking a location on the built-in map.

**⚠️ File name format required:** `LatXX.XXXX_LongYY.YYYY_Recorder_Date_Time.wav`

## Requirements

```r
install.packages(c("shiny", "shinyFiles", "leaflet", "bslib"))
```

## Usage

```r
shiny::runApp("RENAME/0_App_0_Rename.R")
```
or open `0_App_0_Rename.R` in RStudio and click "Run App".

1. Click **"Choose Folder"** to select the directory containing the files to rename.
2. Get a prefix either by clicking on the map (auto lat/long) or by typing it manually — respect the required syntax above.
3. Click **"Preview Changes"** to check the new file names before applying anything.
4. Click **"Apply Renaming"**.

**⚠️ Warning:** renaming is irreversible — original file names cannot be restored afterwards.


