#Load Open Air Trajetories

read.tdumps <- function(hours = 72, pathtraj, height=NULL, location=NULL,
                        pattern,pathtrajout="",onlyseason=NULL,outfilename="trajcombined.txt", 
                        numheaderlines=13) {
  library(dplyr)
  ## find tdump files
  #files <- Sys.glob(paste0(pathtraj,"tdump*"))
  files <-list.files(path = pathtraj, recursive = T, pattern = pattern,full.names = TRUE)
  # if (!is.null(height)==TRUE){ 
  #     files <-grep(paste0("_",height),paste0(pathtraj,files),value=TRUE)}
  if (!is.null(location)==TRUE){ 
      files <-grep(paste0(location),paste0(files),value=TRUE)}
  if (!is.null(onlyseason)==TRUE){ 
      datepattern = paste0("([0-9]{2})([0-9]{4})",location)
      month<-as.numeric(str_match(files,datepattern)[,2])
      seasonnames=c('V', 'O', 'I','P')
      indx <- setNames( rep(seasonnames,each=3), c(12,1:11))
      season <- unname(indx[as.character(month)])
      filesperiod<-subset(files,season==onlyseason)
  } else 
      filesperiod<-files
  
  #files <- as.data.frame(filename)
  output <- file(paste0(pathtrajout,outfilename), 'w')
  for (i in filesperiod){
    # print(i)
    input <- readLines(i)
    numheaderlines=grep("     1 PRESSURE", input)
    input <- input[-c(1:numheaderlines)] # delete header
    writeLines(input, output)
  }
  close(output)
  
  ## read the combined txt file
  traj <- read.table(paste0(pathtrajout, outfilename), header = FALSE)
  traj <- subset(traj, select = -c(V2, V7, V8))
  traj <- rename(traj, c(receptor=V1, year=V3 , month=V4, day=V5,
                         hour=V6, hour.inc=V9, lat=V10 , lon=V11,
                         height=V12, pressure=V13))
  ## hysplit uses 2-digit years ...
  year <- traj$year[1]
  if (year < 50) traj$year <- traj$year + 2000 else traj$year <- traj$year + 1900
  traj$date2 <- with(traj, ISOdatetime(year, month, day, hour, min = 0, sec = 0,
                                       tz = "GMT"))
  ## arrival time
  traj$date <- traj$date2 - 3600 * traj$hour.inc
  traj$datemerge <-as.POSIXct(traj$date , format = "%d%m%Y", tz = "GMT")
  traj
}

add.chem <-function(pollutant, traj){
  traj$datemerge <- traj$date
  traj$datemerge <- format(traj$datemerge,format="%d%m%Y")
  traj$datemerge <- as.POSIXct(strptime(traj$datemerge, format = "%d%m%Y", tz = "GMT"))
  trajconchem<-merge(traj,pollutant)
  trajconchem
}