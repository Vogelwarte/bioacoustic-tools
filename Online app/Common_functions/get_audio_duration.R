

# Function to obtain duration of acoustic files (ffmpeg needs to be installed on the computer for it to work)
get_audio_duration <- function(filepath) {
  # ffprobe command : obtain file duration in seconds (ffmpeg needs to be installed on the computer)
  cmd <- sprintf('ffprobe -v error -show_entries format=duration -of default=noprint_wrappers=1:nokey=1 "%s"', filepath)
  duration <- system(cmd, intern = TRUE)
  # convert to numeric
  as.numeric(duration)
}

