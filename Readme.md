# Bioacoustic tools

This repository hosts Shiny web applications for analyzing, visualizing or optimizing bird detection results from the BirdNET algorithm. These tools provide interactive summaries, spectrograms, audio playback, phenological analysis of bird vocalizations and recording schedule evaluation.

**Note**: Theses applications are independent third-party tools and are not affiliated with or endorsed by the BirdNET developers. They are designed to work with BirdNET output files but are developed and maintained separately.

## Features

- **RENAMING**: Interactive renaming of audio or text files
- (**Data Visualization**: Interactive summary plots showing detection confidence and species distribution)
- **FILTERING**: Interactive summary of birdNET outputs and fitlering of the data
- **SPECTRO**: View spectrograms of detected bird calls with customizable parameters
(- **Audio Playback**: Built-in audio player with simple filtering and noise reduction capabilities)
(- **Data Export**: Export filtered datasets and audio segments in CSV, WAV, or MP3 formats)
- **PHENOLOGY**: Generate temporal activity patterns with sunrise/sunset overlays
- **SCORE**: Sampling Coverage Optimizer for Recording Effort. Evaluate and optimize your current recording schedule
- **Multi-platform Support**: Works on Windows, macOS, and Linux

## Requirements

- R (version 4.1 or higher)
- RStudio (recommended)
- Required R packages (automatically installed via pacman): # Not up to date
  - shiny, bslib, data.table, ggplot2, DT, seewave, tuneR, av
  - tidyr, stringr, stringi, lubridate, hms, parallel
  - mgcv, gratia, suncalc, warbleR, audio, base64enc
- External dependencies:
  - FFmpeg (for MP3/FLAC file support)
  - Audacity (optional, for audio editing)

## Usage

Each App has its own `app.R` file. To run an App, you can use the following command in R: 
```r
shiny::runApp("path/to/BirdNET_ResChecker_app")
```
or open the `app.R` file in RStudio and click "Run App".

### Available Applications

**App_0_Rename:**
A Shiny web application for easily renaming files by adding prefix using the synthax used with these shiny apps (Lat_Long_RECORDER_DATE_TIME).

**App_1_Filtering:**
A Shiny web application for filtering and analyzing BirdNET output files, in relation with species-specific confidence threshold file input. This tool helps researchers and ecologists filter species detections by confidence thresholds, visualize results, and generate summaries with flexible date and GPS options. The App uses combined selection tables from BirdNET.

**App_2_Spectro:**
TO UPDATE: A Shiny web application for analyzing and visualizing bird detection results from the BirdNET algorithm. This tool provides interactive summaries, spectrograms and audio playback of bird vocalizations. Detections validation will be implemented in future versions.

**App_2_Phenology:**
A Shiny web application for visualizing bird detection results from the BirdNET algorithm. This tool provides interactive views of the detections over time.

### Apps location
Apps are located in each folder containing their names. 

### Common Functions
The necessary functions are located in the `R/Common_functions` folder. They will be sourced automatically when running the Apps, as long as the `Common_functions` folder is in the same location as the Shiny `app.R` file.

## Input Files Structure

Expected input files:
- Detection results: `*.BirdNET.selection.table.txt` or compiled selection tables: `*BirdNET_SelectionTable.txt`
- Audio files: `*.wav`, `*.WAV`, `*.mp3`, `*.flac` (in matching folder structure). Integration of other file formats and improvements over the reading time of audio files are in progress.

## Contributors

- Christophe Sahli
- Amandine Serrurier
- Jean-Nicolas Pradervand

## Acknowledgments

This application is designed to work with output from [BirdNET](https://github.com/kahst/BirdNET-Analyzer), a powerful AI-based bird sound identification system developed by the K. Lisa Yang Center for Conservation Bioacoustics at the Cornell Lab of Ornithology.

