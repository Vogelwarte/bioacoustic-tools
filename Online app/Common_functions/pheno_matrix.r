

# pheno_matrix functions code from Amandine Serrurier & Jean-Nicolas Pradervand

# function
pheno_matrix<-function(Voc, SP = "All species", Unit, Confidence1, sunrise, LAT, LONG, UTC, xlim_plot = NULL, nocturnal = F){
  
  # Voc <- DT
  # SP <- "Chouette hulotte"
  # Unit <- 15
  # Confidence1 = 0.8
  # sunrise = T
  # LAT = 46.9
  # LONG = 6.7
  # UTC = "Europe/Paris"
  
    
  #Voc=data_ac
  if(!"All species"%in%SP) {
    Voc <- Voc[Voc$Common.Name%in%SP,] # CHS, species choice (can be a single species or a vector of species)
    
  }
  
  VocRedConf<-Voc[Voc$Confidence>=Confidence1,]#filter by Conf score
  
  #print(paste("nrow(VocRedConf) = ", nrow(VocRedConf)))
  
  if(nrow(VocRedConf)==0){
  return("Error: species not present in the dataset or confidence threshold set too high for the species selected")
  } else {
    
    # title for plots
    if(length(SP)>1){
      SP_title <- paste0("Multiple detected species considered, n = ", length(unique(Voc$Common.Name)), " species")
    } else {
      SP_title <- SP
      if(SP=="All species") {
        SP_title <- paste0("All detected species considered, n = ", length(unique(Voc$Common.Name)), " species")
      }
    }
    
    SeqMin<-seq(0, (1440-Unit), by =Unit)#sequence of minute per time Unit
    
    
    #attr(VocRedConf$Start_segment, "tzone") <- UTC   # = UTC+01:00
    #VocRedConf$Start_segment <- VocRedConf$Start_segment  + hours(UTC)   # = UTC+01:00
    
    # attribute correct UTC
    VocRedConf$time_hms<-as_hms(VocRedConf$Start_segment) # CHS addition: as_hms to solve any format issue
    #attr(VocRedConf$Start_segment, "tzone") <- UTC   # = UTC+01:00
    #VocRedConf$Start_segment <- VocRedConf$Start_segment  + hours(UTC)   # = UTC+01:00
    
    SeqMin_hours<-as_hms(format(as.POSIXct("0:00:00", format="%H:%M:%S")+minutes(SeqMin),format="%H:%M:%S")) #sequence of time per unit # CHS addition: as_hms to solve any format issue
    year1=unique(year(VocRedConf$Start_segment))[1] # CHS modified the datetime source and added [1] to only select a single year. Possible issues: several years and leap year (année bissextile) ?
    SeqDate=data.frame(seq.Date(as.Date(paste0(year1,"-01-01", sep="")), as.Date(paste0(year1,"-12-31", sep="")), by="day")) #sequence of date
    colnames(SeqDate)="date"
    SeqDate$DOY=yday(as.Date(SeqDate$date))#transform date to DOY
    VocRedConf$DOY=yday(as.Date(VocRedConf$Start_segment))#transform date to DOY # CHS modified datetime source
    
    VocMatrix<-as.data.frame(matrix(nrow=length(SeqMin), ncol=nrow(SeqDate)))#empty matrix for number of vocalization per time unit
    colnames(VocMatrix)=seq(1,nrow(SeqDate),by=1)
    
    # CHS Modification of the matrix calculations (vectorization to gain efficiency)
    
    dt <- as.data.table(VocRedConf)# package data.table is very efficient to aggregate date
    # Convert time_hms in seconds ("cut" function hereafter needs numeric values and not hms)
    time_numeric <- as.numeric(dt$time_hms)

    # Convert SeqMin_hours in seconds
    SeqMin_hours_numeric <- as.numeric(SeqMin_hours)
    SeqMin_hours_numeric <- unique(c(SeqMin_hours_numeric, max(time_numeric))) # addition to include all time breaks
    # create a cut according to time_hms / SeqMin_hours
    dt$interval <- cut(time_numeric, breaks = SeqMin_hours_numeric, right = F,include.lowest = T, labels = FALSE)
    # nb of Voc per time Unit
    counts <- dt[, .N, by = .(interval, DOY)]
    
    
    # set base matrix
    VocMatrix <- matrix(0, nrow = length(SeqMin), ncol = nrow(SeqDate),
                        dimnames = list(NULL, SeqDate$DOY))
    
    # Fill matrix
    for (row in 1:nrow(counts)) {
      i <- counts$interval[row]
      j <- which(SeqDate$DOY == counts$DOY[row])
      VocMatrix[i, j] <- counts$N[row]
    }
    
    #transform matrix in DF 
    VocMatrix=as.data.frame(VocMatrix)
    seqY=seq(0, 24-(Unit/60), by=Unit/60)#time proportion 
    matrix=cbind(seqY, VocMatrix)
    matrix=data.matrix(matrix)
    matrix2=matrix%>%as_tibble() %>%gather(key="DOY", value="Vocs", -1)#transform matrix as long data
    matrix2$date_vocmatrix=as.Date(as.numeric(matrix2$DOY), origin = paste0(as.numeric(year1)-1,"-12-31", sep=""))#date in DOY
    midnight<-as.POSIXct("00:00:00", format="%H:%M:%S")
    matrix2$time=format(midnight + (3600*matrix2$seqY), format="%H:%M")#add time per unit
    if (!is.null(xlim_plot)) {
      SeqDateDOY <- SeqDate[SeqDate$date%in%seq(xlim_plot[1], xlim_plot[2], by = "days"),]$DOY # CHS dealing with xlim_plot if chosen in the function arguments
      matrix2<-matrix2[matrix2$DOY%in%SeqDateDOY,]
    }
    print("tibble done")
    
    if (nocturnal==F) {
      a=ggplot()
      
      if (sunrise==TRUE){
        #build a table with sunrise dawn dusk sunset...
        Sun<-getSunlightTimes(date=matrix2$date_vocmatrix, lat=LAT, lon=LONG, keep=c("sunrise", "sunset", "dawn", "dusk"), tz=UTC) #CHS changed SeqDate for matrix2
        #add the right UTC 
        Sun=data.frame(cbind(matrix2$DOY, format(as.POSIXct(Sun$sunset), format="%H:%M:%S"),format(as.POSIXct(Sun$sunrise), format="%H:%M:%S")), 
                       format(as.POSIXct(Sun$dawn), format="%H:%M:%S"),format(as.POSIXct(Sun$dusk), format="%H:%M:%S"))
        Sun$date=matrix2$date_vocmatrix
        colnames(Sun)=c("DOY", "sunset", "sunrise","dawn", "dusk","date")
        noon=format(as.POSIXct("12:00", format="%H:%M"),format="%H:%M")
        midnight2=format(as.POSIXct("23:59",format="%H:%M"),format="%H:%M")
        midnight3=format(as.POSIXct("00:00",format="%H:%M"),format="%H:%M")
        
        #build graph with sunrise and sunset dawn dusk
        a=a+
          geom_ribbon(data=Sun,aes(x=date, ymax=as.POSIXct(format(sunrise), format="%H:%M"), ymin=as.POSIXct(format(midnight3), format="%H:%M")), fill="#CD853F", alpha=0.5)+
          geom_ribbon(data=Sun,aes(x=date, ymin=as.POSIXct(format(sunset), format="%H:%M"), ymax=as.POSIXct(format(midnight2), format="%H:%M")), fill="#CD853F", alpha=0.5)+
          geom_ribbon(data=Sun,aes(x=date, ymax=as.POSIXct(format(dawn), format="%H:%M"), ymin=as.POSIXct(format(midnight3), format="%H:%M")), fill="dodgerblue4", alpha=0.5)+
          geom_ribbon(data=Sun,aes(x=date, ymin=as.POSIXct(format(dusk), format="%H:%M"), ymax=as.POSIXct(format(midnight2), format="%H:%M")), fill="dodgerblue4", alpha=0.5)+
          #geom_ribbon(data=test,aes(x=date, ymax=as.POSIXct(format(sunrise), format="%H:%M"), ymax=as.POSIXct(format(night), format="%H:%M")), fill="white", alpha=1)
          theme(panel.background = element_rect(fill = "white"),panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.background = element_rect(fill = "white"),
                text=element_text(size=18))+
          #geom_tile(matrix2, mapping=aes(date_vocmatrix,as.POSIXct(format(time), format="%H:%M"), fill=Vocs),colour="grey", alpha=0.6)+
          geom_tile(matrix2, mapping=aes(date_vocmatrix,as.POSIXct(format(time), format="%H:%M"), fill=Vocs), color = "grey")+ # plots grid
          geom_tile(matrix2, mapping=aes(date_vocmatrix,as.POSIXct(format(time), format="%H:%M"), fill=Vocs), color = NA)+ # plots no grid so allows gradient to then appear full opacity on the next layer
          scale_fill_gradientn(colours = c("#D9D9D900","#FCFFA4FF","#F98C0AFF" ,"#BB3754FF","#56106EFF","#000004FF"), na.value="#D9D9D900")+
          #scale_x_date(date_labels = "%d%b",date_breaks = "1 month")+
          labs(fill = paste0( "Number of\ndetections/" , Unit, "mins", sep=" "), y="Time", x="Date")+
          scale_y_datetime(date_breaks = "2 hour", date_labels ="%H:%M", expand=c(0,0))+
          ggtitle(label = SP_title)
      }
      else {
        #graph without sunrise and sunset
        a=a+
          geom_tile(matrix2, mapping=aes(date_vocmatrix,as.POSIXct(format(time), format="%H:%M"), fill=Vocs),colour="grey")+ # plots grid
          geom_tile(matrix2, mapping=aes(date_vocmatrix,as.POSIXct(format(time), format="%H:%M"), fill=Vocs),colour=NA)+ # plots no grid so allows gradient to then appear full opacity on the next layer
          scale_fill_gradientn(colours = c("#D9D9D900","#FCFFA4FF","#F98C0AFF" ,"#BB3754FF","#56106EFF","#000004FF"), na.value="#D9D9D900")+
          #scale_x_date(date_labels = "%d%b",date_breaks = "1 month")+
          labs(fill = paste0( "Number of\ndetections/" , Unit, "mins", sep=" "), y="Time", x="Date")+
          scale_y_datetime(date_breaks = "2 hour", date_labels ="%H:%M")+
          theme(text=element_text(size=18))+
          ggtitle(label = SP_title)
        
      }
      
      return(a)
    } else {
      
      # get sunrise and dawn time 
      Sun<-getSunlightTimes(date=matrix2$date_vocmatrix, lat=LAT, lon=LONG, keep=c("sunrise", "sunset", "dawn", "dusk"), tz=UTC) #CHS changed datetime source
      #add the right UTC 
      Sun=data.frame(cbind(matrix2$DOY, format(as.POSIXct(Sun$sunset), format="%H:%M:%S"),format(as.POSIXct(Sun$sunrise), format="%H:%M:%S")), 
                     format(as.POSIXct(Sun$dawn), format="%H:%M:%S"),format(as.POSIXct(Sun$dusk), format="%H:%M:%S"))
      Sun$date=matrix2$date_vocmatrix
      colnames(Sun)=c("DOY", "sunset", "sunrise","dawn", "dusk","date")
      
      
      #add indicator of time for graph 
      today=Sys.Date()#today date for graph
      today2=today+hours(24)#tomorrow date for graph 
      noon=format(as.POSIXct("12:00", format="%H:%M"),format="%H:%M")
      midnight=format(as.POSIXct("23:59",format="%H:%M"),format="%H:%M")
      
      matrix2$graph_date[format(as.POSIXct(matrix2$time, format="%H:%M"),format="%H:%M")>=noon&
                           format(as.POSIXct(matrix2$time, format="%H:%M"),format="%H:%M")<midnight]=format(as.POSIXct(today), format="%Y-%m-%d")
      matrix2$graph_date[format(as.POSIXct(matrix2$time, format="%H:%M"),format="%H:%M")<noon&
                           format(as.POSIXct(matrix2$time, format="%H:%M"),format="%H:%M")<=midnight]=format(as.POSIXct(today)+hours(24), format="%Y-%m-%d")
      
      matrix2$test_graph=as.POSIXct(paste(matrix2$graph_date, matrix2$time), format="%Y-%m-%d %H:%M")
      
      
      
      limits1=as.POSIXct(paste(today, noon), format="%Y-%m-%d %H:%M")
      limits2=as.POSIXct(paste(today2, noon), format="%Y-%m-%d %H:%M")
      
      Sun$graph_date_sunrise=as.POSIXct(paste(today2, Sun$sunrise), format="%Y-%m-%d %H:%M") 
      Sun$graph_date_sunset=as.POSIXct(paste(today, Sun$sunset), format="%Y-%m-%d %H:%M") 
      
      Sun$graph_date_dawn=as.POSIXct(paste(today2, Sun$dawn), format="%Y-%m-%d %H:%M") 
      Sun$graph_date_dusk=as.POSIXct(paste(today, Sun$dusk), format="%Y-%m-%d %H:%M") 
      
      graph_date_midnight=as.POSIXct(paste(today2, "00:00:00"), format="%Y-%m-%d %H:%M")
      graph_date_noon=as.POSIXct(paste(today, "12:00:00"), format="%Y-%m-%d %H:%M") 
      
      # CHS adjustment to fix the temporal slide at midnight occurring otherwise on the plot
      matrix2[format(as.POSIXct(matrix2$time, format="%H:%M"),format="%H:%M")<noon,]$date_vocmatrix <- matrix2[format(as.POSIXct(matrix2$time, format="%H:%M"),format="%H:%M")<noon,]$date_vocmatrix-1
      
      
      a=ggplot()
      if (sunrise==TRUE){
        
        a=a+
          geom_ribbon(Sun, mapping=aes(x=date, ymax=as.POSIXct(graph_date_sunrise,format="%Y-%m-%d %H:%M"), ymin=as.POSIXct(graph_date_midnight, format="%Y-%m-%d %H:%M")), fill="#CD853F", alpha=0.5)+
          geom_ribbon(Sun, mapping=aes(x=date, ymin=as.POSIXct(graph_date_sunset,format="%Y-%m-%d %H:%M"), ymax=as.POSIXct(graph_date_midnight, format="%Y-%m-%d %H:%M")), fill="#CD853F", alpha=0.5)+
          geom_ribbon(Sun, mapping=aes(x=date, ymax=as.POSIXct(graph_date_dawn,format="%Y-%m-%d %H:%M"), ymin=as.POSIXct(graph_date_midnight, format="%Y-%m-%d %H:%M")), fill="dodgerblue4", alpha = 0.5)+
          geom_ribbon(Sun, mapping=aes(x=date, ymin=as.POSIXct(graph_date_dusk,format="%Y-%m-%d %H:%M"), ymax=as.POSIXct(graph_date_midnight, format="%Y-%m-%d %H:%M")), fill="dodgerblue4", alpha = 0.5)+
          geom_tile(matrix2, mapping=aes(date_vocmatrix,as.POSIXct(format(test_graph), format="%Y-%m-%d %H:%M"), fill=Vocs),color="grey")+ # plots grid
          geom_tile(matrix2, mapping=aes(date_vocmatrix,as.POSIXct(format(test_graph), format="%Y-%m-%d %H:%M"), fill=Vocs),color=NA)+ # plots no grid so allows gradient to then appear full opacity on the next layer
          scale_fill_gradientn(colours = c("#D9D9D900","#FCFFA4FF","#F98C0AFF" ,"#BB3754FF","#56106EFF","#000004FF"), na.value="#D9D9D900")+
          scale_y_datetime("Time", #change axis title back to time
                           #Set limits to show a full day
                           limits = c(limits1, limits2),
                           expand = c(0,0), #turn off axis limit expansion
                           date_breaks = "2 hours", #set custom axis breaks
                           date_labels = "%H:%M")+
          theme(panel.background = element_rect(fill = "white"),panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.background = element_rect(fill = "white"),
                text=element_text(size=18))+
          #scale_x_date(date_labels = "%d%b",date_breaks = "1 month")+
          labs(fill = paste0( "Number of\ndetections/" , Unit, "mins", sep=" "), y="Time", x="Date")+
          ggtitle(label = SP_title)
        
      }
      else {
        a=a+
          geom_tile(matrix2, mapping=aes(date_vocmatrix,as.POSIXct(format(test_graph), format="%Y-%m-%d %H:%M"), fill=Vocs),color="grey")+ # plots grid
          geom_tile(matrix2, mapping=aes(date_vocmatrix,as.POSIXct(format(test_graph), format="%Y-%m-%d %H:%M"), fill=Vocs),color = NA)+ # plots no grid so allows gradient to then appear full opacity on the next layer
          scale_fill_gradientn(colours = c("#D9D9D900","#FCFFA4FF","#F98C0AFF" ,"#BB3754FF","#56106EFF","#000004FF"), na.value="#D9D9D900")+
          scale_y_datetime("Time", #change axis title back to time
                           #Set limits to show a full day
                           limits = c(limits1, limits2),
                           expand = c(0,0), #turn off axis limit expansion
                           date_breaks = "2 hours", #set custom axis breaks
                           date_labels = "%H:%M")+
          theme(panel.background = element_rect(fill = "white"),panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.background = element_rect(fill = "white"),
                text=element_text(size=18))+
          #scale_x_date(date_labels = "%d%b",date_breaks = "1 month")+
          labs(fill = paste0( "Number of\ndetections/" , Unit, "mins", sep=" "), y="Time", x="Date")+
          ggtitle(label = SP_title)
      }
      
      return(a)
    }
    
  }
  
}


