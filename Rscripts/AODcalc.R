library(readr)
library(readxl)
library(openair)
library(ggplot2)
library(dplyr)
setwd("~/Documents/paper_facu/Rscripts")
CEILAP_BA_hora <- read_csv("../data/DatosCEILAP/20190101_20200630_CEILAP-BA.lev20", 
                      col_types = cols(`Date(dd:mm:yyyy)` = col_datetime(format = "%d:%m:%Y")), 
                      skip = 6, na = "-999.000000")
CEILAP_BA_hora$date = as.POSIXct(paste(CEILAP_BA_hora$`Date(dd:mm:yyyy)`, CEILAP_BA_hora$`Time(hh:mm:ss)`), format="%Y-%m-%d %H:%M:%S",tz='UTC')
CEILAP_BA_hora$date <- as.POSIXct(CEILAP_BA_hora$date, tz='UTC')
CEILAP_BA_hora<-timeAverage(CEILAP_BA_hora,avg.time = "hour")
CEILAP_BA_hora$filterday <- as.Date(CEILAP_BA_hora$date - 12*3600)

BA_events_testM <- read_excel("../BA_events_testMnew.xlsx")
BA_events_testM$date <- as.POSIXct(BA_events_testM$date, tz='UTC')
BA_events_testM$filterday <- as.Date(BA_events_testM$date)
BA_events_testM$Event_M <- as.factor(BA_events_testM$Event_M)
BA_events_testM$Event_F <- as.factor(BA_events_testM$Event_F)
CEILAP_BA_hora2= merge(CEILAP_BA_hora, BA_events_testM, by="filterday", suffixes = c("_CEILAPBA","_sample"))


# CEILAP_BA_hora2 = merge(CEILAP_BA_hora, BA_events_testM, by.x="Date(dd:mm:yyyy)", by.y="date")
CEILAP_BA_hora3=CEILAP_BA_hora2[,c("date","AOD_440nm","440-870_Angstrom_Exponent","ID","Event_F","OrigenTag")]

png("../images/AOD_FIREOrigins_colors_vf.png", res=300, height=1.5*800,width = 1.5*1000)
ggplot()+  geom_point(data = CEILAP_BA_hora2, 
                     aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`), color="grey50",fill="grey50",shape=21, size=2)+
  geom_point(data = CEILAP_BA_hora2[(CEILAP_BA_hora2$Event_F %in% c("SI","SO")) & (CEILAP_BA_hora2$AOD_440nm>=0.16) ,], 
             aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`,color=OrigenTag,fill=OrigenTag), shape=21, size=2)+
  scale_color_discrete(name = "Event Origin",na.translate = F)+scale_fill_discrete(name = "Event Origin",na.translate = F)+
  labs(y=expression(alpha["440-880nm"]), x=expression(AOD["440nm"]))+
  theme_bw()+ theme(legend.position="top")+ geom_vline(xintercept=0.1613,color="black",linetype="dashed")+
  guides(color = guide_legend(nrow = 2)) + 
  scale_x_continuous(
    breaks = function(x) unique(c(0.1613, scales::pretty_breaks()(x))),
    labels = function(x) {
      ifelse(x == 0.1613, paste0("\n", format(x, nsmall = 2, digits = 2)), format(x, nsmall = 2, digits = 2))
    }
  )
dev.off()
ggplot()+ geom_point(data = CEILAP_BA_hora2, 
                     aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`))+
  geom_point(data = CEILAP_BA_hora2[CEILAP_BA_hora2$OrigenTag %in% c("AUS"),], 
             aes(x=AOD_440nm,y=`440-870_Angstrom_Exponent`,color=OrigenTag))+
  geom_hline(yintercept=1.13,color="grey50")+scale_color_discrete(name = "Event Origin")+
  geom_vline(xintercept=0.1613,color="grey50")+
  theme_bw()
BOLCEILAP_BA_hora2$Day_of_Year=factor(CEILAP_BA_hora2$Day_of_Year)
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

CEILAP_BA <- read_csv("../data/DatosCEILAP/20190101_20200630_CEILAP-BA.lev20", 
                      col_types = cols(`Date(dd:mm:yyyy)` = col_datetime(format = "%d:%m:%Y")), 
                      skip = 6, na ="-999.000000")

CEILAP_BA=CEILAP_BA[,c(1,2,22,65)]

CEILAP_BA$date = as.POSIXct(paste(CEILAP_BA$`Date(dd:mm:yyyy)`, CEILAP_BA$`Time(hh:mm:ss)`), format="%Y-%m-%d %H:%M:%S",tz='UTC')
CEILAP_BA$date <- as.POSIXct(CEILAP_BA$date, tz='UTC')
CEILAP_BA$date = CEILAP_BA$date - 12*3600
CEILAP_BA=timeAverage(CEILAP_BA,avg.time = "day",na.rm=T)
CEILAP_BA = CEILAP_BA[,c(1,4,5)]