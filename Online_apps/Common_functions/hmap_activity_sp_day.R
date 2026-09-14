
##### HEATMAP ACTIVITY SPECIES AND DAY #####
# create by Amandine Serrurier with recording from PVPs projects Vogelwarte - updated 29.09.2025

# the code will count detections of species (filtered by confidence score or not) by day for each species
#it will then create a heatmap of activity along the season for each species 
# Original code: ~Original scripts/heatmap_activity_SP&Day_ASe.R

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


##### split columns #####
#depending on your folder organisation, you have to define which "parts" correspond to the information you need:
# my usual path is " "/Volumes/PortableSSD/Grengiols/GRE_A/3/Data/SMA03227_20240726_144657.GPU.BirdNET.selection.table.txt"
#which means my site name is always placed 4th position, my recorder name 5th position
#file name of recording is always the basename

parts <- str_split(big_table$FileNameCC, "/", simplify = TRUE) #split parts of your path to extract information
pos_SITE= 4
pos_REC= 5

format_data<- big_table %>%
  mutate(
    file = basename(FileNameCC),
    recorder = parts[,pos_REC], #recorder name
    site     = parts[,pos_SITE], #site name
    datetime_str = str_extract(file, "\\d{8}_\\d{6}"), #datetime recorder format
    datetime = ymd_hms(datetime_str, tz = "UTC"), #datetime UTC format
    date = as.Date(datetime)
  )

##### format data for graph #####
#Count your detection by date and by species
species_date=format_data%>% group_by(date, Name) %>% summarise(Count = n())

#Standardise your count between 0 and 1 per species:
df_norm <- species_date %>%
  group_by(Name) %>%
  mutate(Activite = (Count - min(Count)) / (max(Count) - min(Count))) %>%
  ungroup() %>% mutate(Name = factor(Name, levels = rev(sort(unique(Name)))))

#add an option to filter by species, by confidence score and by date would be nice !

##### VISUALISATION plot the result #####
ggplot(df_norm, aes(x=as.Date(date), y=Name, fill=Activite))+geom_tile(color="grey", alpha=0.9)+
  scale_fill_gradient(low = "#E0F2F1", high = "#00695C", na.value="#00695C" )+
  geom_hline(yintercept = seq(0.5, length(unique(df_norm$Name)) + 0.5, 1), 
             color = "grey", alpha = 0.9, linewidth=0.3)+
  theme(panel.background = element_rect(fill = "white",),panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor =element_blank(),
        plot.background = element_rect(fill = "white"))+
  labs(fill = "Activity\nstandardised\n", y="Species", x="Date")+scale_x_date(date_breaks = "1 weeks", date_labels = "%d %b")




















