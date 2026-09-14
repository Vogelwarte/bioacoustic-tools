
##### HEATMAP RECORDING TIME #####
# create by Amandine Serrurier with recording from PVPs projects Vogelwarte - updated 29.09.2025

# the code will list of the recording available for one or several recorder, identify time and date avaialble, 
# and create a graph along the season with recording time. Allows for exemple to identify missing data 
# Original file: ~Original scripts/heatmap_recording_ASe.R


##### packages ##### must be removed when function works for the app
library(dplyr)
library(stringr)
library(lubridate)
library(ggplot2)
library(readr)      
library(tidyr)      
library(suncalc)
library(splusTimeDate)
library(fs)

##### load data #####
rec_path="/Users/amandineserrurier/Documents/VoWa/PVPs" #path with recordings
file_list=dir_ls(rec_path, regexp = "wav", recurse = TRUE, fail=FALSE) #list filenames of recording
path_table=data.frame(file) #create a table with filepaths
path_table =  path_table %>% rename(path = 1)  # name path column of table "path"

#df_25=read.csv2("/Users/amandineserrurier/Downloads/summary_file_250925.csv")

##### extract datetime #####
parts <- str_split(path_table$path, "/", simplify = TRUE) #split parts of your path to extract information

#depending on your folder organisation, you have to define which "parts" correspond to the information you need:
# my usual path is "Z:/SciData/DOM_Forschung/PAM/transfer/audio_recordings/460410_Photovoltaikanlagen/2025/Belalp/BEL_A/1/Data/SMA03365_20250722_085803.wav"
#which means my site name is always placed 9th position, my recorder name 10th position
#file name of recording is always the basename

pos_SITE= 9
pos_REC= 10

path_table <- path_table %>%
  mutate(
    file = basename(path),
    recorder = parts[,pos_REC], #recorder name
    site     = parts[,pos_SITE], #site name
    datetime_str = str_extract(file, "\\d{8}_\\d{6}"), #datetime recorder format
    datetime = ymd_hms(datetime_str, tz = "UTC"), #datetime UTC format
    date = as.Date(datetime)
  )

# you can choose here the recorder or site you want to display or juste display everything
pkSITE="Belalp"
#or 
pkREC= "GIB_A"

path_table=path_table[path_table$site == pkSITE,]

##### summarise recording time #####
rec_day<- path_table%>%
  group_by(recorder, date, site ) %>%
  summarise(hours = n(), .groups = "drop")

# Define all dates, recorders and site for graph
all_dates <- seq(min(na.omit(rec_day$date)), max(na.omit(rec_day$date)), by = "day") #all date
recorders<- unique(rec_day$recorder) #list recorders
sites<- unique(rec_day$site) #list sites

# fill long table with all possibles dates and recorders
df_complete <- expand.grid(recorder = recorders, date = all_dates, site=sites) %>%
  left_join(rec_day, by = c("recorder", "date") ) %>%
  mutate(hours = replace_na(hours, 0))%>% rename(location="site.x") %>%select(-site.y)

##### VISUALISATION tile graph #####
# 6. Graph heatmap
# df_complete$hours <- ifelse(df_complete$hours == 0,  NA, df_complete$hours)
ggplot(df_complete) +
  geom_tile(aes(x = date, y = recorder, fill = hours), color = "grey90", linewidth = 0.05, linejoin = "round") +
  # geom_tile(aes(x = as.Date(rec_end_25), y = recorder),color = "grey90",  fill="blue", linewidth = 0.05) +
  # geom_tile(aes(x = as.Date(rec_deployed_2025), y = recorder), fill="blue" , color = "grey90", linewidth = 0.05) +
  scale_fill_continuous(low="white", high="darkred", na.value="transparent") +
  theme_minimal() +
  labs(
    x = "Date",
    y = "Recorders",
    title = "Number of recording per day and recorder"
  ) +
  
  scale_x_date(date_breaks = "2 days", date_labels = "%b %d",expand = expansion(add = 10) )+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.y = element_blank(),  panel.grid.major.x  = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())































