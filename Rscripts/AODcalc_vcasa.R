library(readr)
library(readxl)
library(openair)
# CEILAP_BA_TOT <- read_csv("data/DatosCEILAP/20190101_20200630_CEILAP-BA.tot_lev20", 
#                           col_types = cols(`Date(dd:mm:yyyy)` = col_datetime(format = "%d:%m:%Y")), 
#                           skip = 6)
CEILAP_BA_hora <- read_csv("data/DatosCEILAP/20190101_20200630_CEILAP-BA.lev20", 
                      col_types = cols(`Date(dd:mm:yyyy)` = col_datetime(format = "%d:%m:%Y")), 
                      skip = 6, na = "-999")
CEILAP_BA <- read_csv("data/DatosCEILAP/20190101_20200630_CEILAP-BA.lev20", 
                      col_types = cols(`Date(dd:mm:yyyy)` = col_datetime(format = "%d:%m:%Y")), 
                      skip = 6, na ="-999.000000")

CEILAP_BA=CEILAP_BA[,c(1,2,22,65)]
CEILAP_BA$date = as.POSIXct(paste(CEILAP_BA$`Date(dd:mm:yyyy)`, CEILAP_BA$`Time(hh:mm:ss)`), format="%Y-%m-%d %H:%M:%S",tz='UTC')
CEILAP_BA$date = CEILAP_BA$date - 12*3600
subsets=subset(CEILAP_BA, CEILAP_BA$AOD_440nm>0.16)
CEILAP_BA=timeAverage(CEILAP_BA,avg.time = "day")
subsetdia=timeAverage(subsets,avg.time = "day",na.rm=TRUE,statistic = "max")
CEILAP_BA = CEILAP_BA[,c(1,4,5)]

BA_events_testM <- read_excel("BA_events_testM.xlsx")
BA_events_testM$date <- as.POSIXct(BA_events_testM$date, tz='UTC')
BA_events_testM$Event_M <- as.factor(BA_events_testM$Event_M)
CEILAP_BA = merge(CEILAP_BA, BA_events_testM, by="date")
CEILAP_BA$diffAOD=CEILAP_BA$AOD_440nm-CEILAP_BA$AOD440
subsetmergeo = merge(subsetdia, BA_events_testM, by="date")

subset1=subset(subsetmergeo,subsetmergeo$AOD440>0.16)
subset2=subset(subset1,subset1$Event_M %in% c("S","SN", "SC","SL"))
subset3=subset(subsetmergeo,subsetmergeo$Event_M %in% c("S","SN", "SL"))
