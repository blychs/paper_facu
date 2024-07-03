library(readr)
library(readxl)
library(openair)
library(ggplot2)
library(dplyr)
# CEILAP_BA_TOT <- read_csv("data/DatosCEILAP/20190101_20200630_CEILAP-BA.tot_lev20", 
#                           col_types = cols(`Date(dd:mm:yyyy)` = col_datetime(format = "%d:%m:%Y")), 
#                           skip = 6)
CEILAP_BA_hora <- read_csv("../data/DatosCEILAP/20190101_20200630_CEILAP-BA.lev20", 
                      col_types = cols(`Date(dd:mm:yyyy)` = col_datetime(format = "%d:%m:%Y")), 
                      skip = 6, na = "-999.000000")

CEILAP_BA <- read_csv("../data/DatosCEILAP/20190101_20200630_CEILAP-BA.lev20", 
                      col_types = cols(`Date(dd:mm:yyyy)` = col_datetime(format = "%d:%m:%Y")), 
                      skip = 6, na ="-999.000000")

CEILAP_BA=CEILAP_BA[,c(1,2,22,65)]

CEILAP_BA$date = as.POSIXct(paste(CEILAP_BA$`Date(dd:mm:yyyy)`, CEILAP_BA$`Time(hh:mm:ss)`), format="%Y-%m-%d %H:%M:%S",tz='UTC')
CEILAP_BA$date <- as.POSIXct(CEILAP_BA$date, tz='UTC')
CEILAP_BA$date = CEILAP_BA$date - 12*3600
CEILAP_BA=timeAverage(CEILAP_BA,avg.time = "day",na.rm=T)
CEILAP_BA = CEILAP_BA[,c(1,4,5)]

BA_events_testM <- read_excel("../BA_events_testMnew.xlsx")
BA_events_testM$date <- as.POSIXct(BA_events_testM$date, tz='UTC')
BA_events_testM$Event_M <- as.factor(BA_events_testM$Event_M)
BA_events_testM$Event_F <- as.factor(BA_events_testM$Event_F)
CEILAP_BA2 = merge(CEILAP_BA, BA_events_testM, by="date")
# CEILAP_BA$diffAOD=CEILAP_BA$AOD_440nm-CEILAP_BA$AOD440

CEILAP_BA_hora2 = merge(CEILAP_BA_hora, BA_events_testM, by.x="Date(dd:mm:yyyy)", by.y="date")

# ggplot()+ geom_point(data = CEILAP_BA2, 
#                      aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`))+
#   geom_point(data = CEILAP_BA2[CEILAP_BA2$Event_M %in% c("S","SL","SN"),], 
#              aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`,color=Event_M))+
#   geom_hline(yintercept=1.13,color="grey50")+
#   geom_vline(xintercept=0.1613,color="grey50")+
#   theme_bw()
# 
# # summary clima
# CEILAP_BA_hist <- read_csv("../data/DatosCEILAP/20100101_20201231_CEILAP-BA.lev20", 
#                            col_types = cols(`Date(dd:mm:yyyy)` = col_datetime(format = "%d:%m:%Y")), 
#                            skip = 6, na = c("-999","-999.000000"))
# summary(CEILAP_BA_hist$AOD_440nm)
# 
# ggplot()+ geom_point(data = CEILAP_BA_hora2, 
#                      aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`))+
#   geom_point(data = CEILAP_BA_hora2[CEILAP_BA_hora2$Event_F %in% c("SI","SO","SN","SP"),], 
#              aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`,color=Event_F))+
#   geom_hline(yintercept=1.13,color="grey50")+
#   geom_vline(xintercept=0.1613,color="grey50")+
#   theme_bw()
# 
# ggplot()+ geom_point(data = CEILAP_BA_hora2, 
#                      aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`))+
#   geom_point(data = CEILAP_BA_hora2[CEILAP_BA_hora2$Event_F %in% c("SI","SO","SN","SP"),], 
#              aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`,color=Local))+
#   geom_hline(yintercept=1.13,color="grey50")+
#   geom_vline(xintercept=0.1613,color="grey50")+
#   theme_bw()

ggplot()+ geom_point(data = CEILAP_BA_hora2, 
                     aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`))+
  geom_point(data = CEILAP_BA_hora2[CEILAP_BA_hora2$Event_F %in% c("SI"),], 
             aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`,color=Tag))+
  geom_hline(yintercept=1.13,color="grey50")+
  geom_vline(xintercept=0.1613,color="grey50")+
  theme_bw()


beCEILAP_BA_hora2$Day_of_Year=factor(CEILAP_BA_hora2$Day_of_Year)
CEILAP_BA_hora2$datefactor=factor(CEILAP_BA_hora2$`Date(dd:mm:yyyy)`)

ggplot()+ geom_point(data = CEILAP_BA_hora2, 
                     aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`))+
  geom_point(data = CEILAP_BA_hora2[CEILAP_BA_hora2$Event_M %in% c("S","SL","SN"),], 
             aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`,color=datefactor))+
  geom_hline(yintercept=1.13,color="grey50")+
  geom_vline(xintercept=0.1613,color="grey50")+
  theme_bw()

ggplot()+ geom_point(data = CEILAP_BA_hora2, 
                     aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`))+
  geom_point(data = CEILAP_BA_hora2[CEILAP_BA_hora2$Event_F %in% c("SI","SF","SN","SP"),], 
             aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`,color=Origen))+
  geom_hline(yintercept=1.13,color="grey50")+
  geom_vline(xintercept=0.1613,color="grey50")+
  theme_bw()

ggplot()+ geom_point(data = CEILAP_BA_hora2, 
                     aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`))+
  geom_point(data = CEILAP_BA_hora2[CEILAP_BA_hora2$Event_F %in% c("SI") ,], 
             aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`,color=Origen))+
  geom_hline(yintercept=1.13,color="grey50")+
  geom_vline(xintercept=0.1613,color="grey50")+
  theme_bw()


ggplot()+ geom_point(data = CEILAP_BA2, 
                     aes(x=`AOD440nm1.5`,y=`Alpha1.5`))+
  geom_point(data = CEILAP_BA2[CEILAP_BA2$Event_F %in% c("SI","SF","SN","SP"),], 
             aes(x=`AOD440nm1.5`,y=`Alpha1.5`,color=Origen))+
  geom_hline(yintercept=1.13,color="grey50")+
  geom_vline(xintercept=0.1613,color="grey50")+
  theme_bw()


ggplot()+ geom_point(data = CEILAP_BA2, 
                     aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`))+
  geom_point(data = CEILAP_BA2[CEILAP_BA2$Event_F %in% c("SI"),], 
             aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`,color=Origen))+
  geom_hline(yintercept=1.13,color="grey50")+
  geom_vline(xintercept=0.1613,color="grey50")+
  theme_bw()

ggplot()+ geom_point(data = CEILAP_BA2, 
                     aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`))+
  geom_point(data = CEILAP_BA2[CEILAP_BA2$Event_F %in% c("SI","SF","SN","SP"),], 
             aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`,color=Origen))+
  geom_hline(yintercept=1.13,color="grey50")+
  geom_vline(xintercept=0.1613,color="grey50")+
  theme_bw()
subsetAOD=subset(CEILAP_BA2,CEILAP_BA2$AOD_440nm>0.16 |CEILAP_BA2$AOD440nm1.5>0.16)
subsetEventos=subset(CEILAP_BA2,CEILAP_BA2$Event_F %in% c("SI"))

                     