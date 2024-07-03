library(openair)
library(ggplot2)
library(lubridate)
library(readxl)
setwd("~/Documents/paper_facu")
pathgraphs="Figures"
upper_windspeed=8.5
colorset="plasma"
# En data_every_hour_obs.csv solo estan actualizados los datos de meteo el resto del archivo es viejo
data <- read.csv("data/data_every_hour_obs_eventM.csv")                           
data$date <- as.POSIXct(data$date, tz='UTC')
mergeddata <- data[,-c(29:40)]
mergeddata$Sb_ng=mergeddata$Sb*1000
mergeddata$As_ng=mergeddata$As*1000
mergeddata$nssK=mergeddata$K-0.6*mergeddata$Fe-0.037*mergeddata$Na
mergeddata$nssK_OC=mergeddata$nssK/mergeddata$C.Orgánico
mergeddata$nssK_EC=mergeddata$nssK/mergeddata$C.Elemental
mergeddata$OC_EC=mergeddata$C.Orgánico/mergeddata$C.Elemental
mergeddata=dplyr::rename(mergeddata,Event= Event_M)
# data$day<-floor_date(data$date-3600*12,"day")
# meteodata <- data[,c(1,47,48,49)]
# # de aca sacar Event_M
# BA_events_testM <- read_excel("BA_events_testM.xlsx", 
#                               col_types = c("date", "numeric", "numeric", 
#                                             "skip", "text"))
# BA_events_testM$date <- as.POSIXct(BA_events_testM$date, tz='UTC')
# mergeddata=merge(meteodata,BA_events_testM, by.x = "day", by.y = "date")
# #De acá saco los valores de PMF
# PMF_BA_full <- read_excel("PMF_BA_full.xlsx", 
#                           sheet = "CONC", col_types = c("date", 
#                                                         "numeric", "numeric", "numeric", 
#                                                         "numeric", "numeric", "numeric", 
#                                                         "skip", "numeric", "numeric", "numeric", 
#                                                         "numeric", "numeric", "numeric", 
#                                                         "numeric", "numeric", "numeric", 
#                                                         "numeric", "numeric", "numeric", 
#                                                         "numeric", "numeric", "numeric", 
#                                                         "numeric", "numeric", "numeric", 
#                                                         "numeric", "numeric", "numeric", 
#                                                         "numeric", "skip", "skip", "skip", 
#                                                         "skip", "skip", "numeric", "skip", 
#                                                         "skip", "skip", "skip", "skip", "skip", 
#                                                         "numeric", "numeric", "numeric"))
# mergeddata=merge(mergeddata,PMF_BA_full,by.x = "day", by.y = "date")
rm(data)
# rm(meteodata)

mergeddata$OC_EC = mergeddata$`C.Orgánico`/mergeddata$`C.Elemental`
mergeddata$OC_K=mergeddata$`C.Orgánico`/mergeddata$K

PPPM25<-polarPlot(mergeddata, pollutant = "PM2.5", statistic = "mean",  min.bin = 2, 
                  upper =upper_windspeed, key.footer="[ug/m3]",key.header = "PM2.5",mis.col = "transparent",
                  cols="inferno")
#polarPlot(mergeddata, pollutant = "OC_EC", statistic = "mean",  min.bin = 2, main='OC/EC')
ppEC<-polarPlot(mergeddata, pollutant = "C.Elemental", statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, key.footer="[ug/m3]",key.header = "EC",mis.col = "transparent",
          cols=colorset)
ppOC<-polarPlot(mergeddata, pollutant = "C.Orgánico", statistic = "mean",  min.bin = 2, 
                upper =upper_windspeed, key.footer="[ug/m3]",key.header = "OC",mis.col = "transparent",
                cols=colorset)
ppOM<-polarPlot(mergeddata, pollutant = "OM", statistic = "mean",  min.bin = 2, 
                upper =upper_windspeed, key.footer="[ug/m3]",key.header = "OM",mis.col = "transparent",
                cols=colorset)
png(paste0(pathgraphs,"/polarplotECOCenmug_",colorset,".png"), width = 717 * 5, height = 339* 5, res = 300)
print(ppEC$plot,split=c(1, 1, 2, 1))
print(ppOC$plot,split=c(2, 1, 2, 1), newpage=FALSE)
dev.off()

png(paste0(pathgraphs,"/polarplotECOCOM_",colorset,".png"), width = 900 * 4, height = 270* 4, res = 300)
print(ppEC$plot,split=c(1, 1, 3, 1))
print(ppOC$plot,split=c(2, 1, 3, 1), newpage=FALSE)
print(ppOM$plot,split=c(3, 1, 3, 1), newpage=FALSE)
dev.off()

ppAs<-polarPlot(mergeddata, pollutant = "As", statistic = "mean",  min.bin = 2, 
                upper =upper_windspeed, key.header="As",key.footer="[ug/m3]",
                mis.col = "transparent", cols=colorset)
ppSb<-polarPlot(mergeddata, pollutant = "Sb", statistic = "mean",  min.bin = 2, 
                upper =upper_windspeed, key.header ="Sb",key.footer="[ug/m3]",
                mis.col = "transparent", cols=colorset)
png(paste0(pathgraphs,"/polarplotAsSb_",colorset,".png"), width = 717 * 5, height = 339* 5, res = 300)
print(ppAs$plot,split=c(1, 1, 2, 1))
print(ppSb$plot,split=c(2, 1, 2, 1), newpage=FALSE)
dev.off()

ppZn<-polarPlot(mergeddata, pollutant = "Zn", statistic = "mean",  min.bin = 2, 
                upper =upper_windspeed, key.header="Zn",key.footer="[ug/m3]",
                mis.col = "transparent", cols=colorset)


ppAs<-polarPlot(mergeddata, pollutant = "As", statistic = "mean",  min.bin = 2, 
                upper =upper_windspeed, key.header="As",key.footer="[ug/m3]",
                mis.col = "transparent")
ppSb<-polarPlot(mergeddata, pollutant = "Sb", statistic = "mean",  min.bin = 2, 
                upper =upper_windspeed, key.header ="Sb",key.footer="[ug/m3]",
                mis.col = "transparent")
png(paste0(pathgraphs,"/polarplotAsSb.png"), width = 717 * 5, height = 339* 5, res = 300)
print(ppAs$plot,split=c(1, 1, 2, 1))
print(ppSb$plot,split=c(2, 1, 2, 1), newpage=FALSE)
dev.off()

ppAs<-polarPlot(mergeddata, pollutant = "K", statistic = "mean",  min.bin = 2, 
                upper =upper_windspeed, key.header="K",key.footer="[ug/m3]",
                mis.col = "transparent", cols=colorset)
ppAs<-polarPlot(mergeddata, pollutant = "As_ng", statistic = "mean",  min.bin = 2, 
                upper =upper_windspeed, key.header="As",key.footer="[ng/m3]",
                mis.col = "transparent", cols=colorset)
ppSb<-polarPlot(mergeddata, pollutant = "Sb_ng", statistic = "mean",  min.bin = 2, 
                upper =upper_windspeed, key.header ="Sb",key.footer="[ng/m3]",
                mis.col = "transparent", cols=colorset)
png(paste0(pathgraphs,"/polarplotAsSbenng_",colorset,".png"), width = 717 * 5, height = 339* 5, res = 300)
print(ppAs$plot,split=c(1, 1, 2, 1))
print(ppSb$plot,split=c(2, 1, 2, 1), newpage=FALSE)
dev.off()


ppNi<-polarPlot(mergeddata, pollutant = "Ni", statistic = "mean",  min.bin = 2, 
                upper =upper_windspeed, key.header="Ni",key.footer="[ug/m3]",
                mis.col = "transparent", cols=colorset)
ppV<-polarPlot(mergeddata, pollutant = "V", statistic = "mean",  min.bin = 2, 
                upper =upper_windspeed, key.header ="V",key.footer="[ug/m3]",
                mis.col = "transparent", cols=colorset)
png(paste0(pathgraphs,"/polarplotNiV_",colorset,".png"), width = 717 * 5, height = 339* 5, res = 300)
print(ppNi$plot,split=c(1, 1, 2, 1))
print(ppV$plot,split=c(2, 1, 2, 1), newpage=FALSE)
dev.off()

# scatterPlot(PMF_BA_full_so,x="NH4",y="nssSO4")
# ## codigo anterior Pablo
# keys = names(data)[2:45]
# print(keys)
# 
# data_period2 = data[which(data$date >= as.POSIXct("2019-05-22", format="%Y-%m-%d") & data$date <= as.POSIXct("2019-06-04", format="%Y-%m-%d")), ]
# print(data_period2)
# windRose(data)#, poll = 'Na.sol', stati = 'cpf', percentile=c(0, 100),  main=paste('Na sol','cpf'))
# 
# #data[date(data$date) > as.POSIXct('2019-05-23') & date(data$date) < as.POSIXct('2019-06-04') ,-1] = NA
# #data[date(data$date) == as.POSIXct(c('2020-02-21', '2020-02-22')),] = NA
# max_perc = 100
# min_perc = 75
# for (title in keys) {
# print(title)
# png(paste('images/median/bv_median_',title,'.png', sep=''))
# polarPlot(data, poll = title, stati = 'median',  main=paste(title,'median'))
# dev.off()
# }
# 
# 
# polarPlot(data)
# 
# 
# data_event = data[which(data$Event%in%c('S', 'SP', 'SN')),]
# max_perc = 25
# min_perc = 0
# for (title in keys) {
#   print(title)
#   png(paste('images/CPF/CPF_Q1_',title,'.png', sep=''))
#   polarPlot(data, poll = title, stati = 'cpf', percentile = c(min_perc,max_perc), main=paste(title,'Q1'))
#   dev.off()
# }
# 
# data_no_event = data[which(!data$Event%in%c('S', 'SP', 'SN')),]
# max_perc = 100
# min_perc = 75
# for (title in keys) {
#   print(title)
#   png(paste('images/CPF/CPF_No_Event_P75_',title,'.png', sep=''))
#   polarPlot(data_no_event, poll = title, stati = 'cpf', percentile = c(min_perc,max_perc), main=paste(title,'P75'))
#   dev.off()
# }
# 
# date <- as.POSIXct("2019-17-7", format="%Y-%m-%d")
# data_period1 = data[which(data$date < as.POSIXct("2019-07-17", format="%Y-%m-%d")),]
# data_period2 = data[which(data$date >= as.POSIXct("2019-07-17", format="%Y-%m-%d") & data$date < as.POSIXct("2019-11-14", format="%Y-%m-%d")), ]
# print(data_period2)
# 
# 
# 
# print("Periods")
# for (title in keys) {
#   print(paste(title, "Period 1"))
#   png(paste("images/CPF/Na_Period1/CPF_Period1_", title, "_Q4.png", sep=""))
#   polarPlot(data_period1, poll = title, stati = 'cpf', percentile = c(75, 100), main=paste(title, 'before Jul17 Q4'))
#   dev.off()
#   print('done Q4')
#   if (title=='Co'){
#     png(paste("images/CPF/Na_Period1/CPF_Period1_", title, "_Q1.png", sep=""))
#     polarPlot(data_period1, poll = title, stati = 'cpf', percentile = c(0, 25), main=paste(title, 'before Jul17 Q1'), limits=c(0, 8.441e-05))
#     dev.off()
#   } else {
#     png(paste("images/CPF/Na_Period1/CPF_Period1_", title, "_Q1.png", sep=""))
#     polarPlot(data_period1, poll = title, stati = 'cpf', percentile = c(0, 25), main=paste(title, 'before Jul17 Q1'))
#     dev.off()
#   }
# 
#   print('done Q1')
#   print(paste(title, "Period 2"))
#   png(paste("images/CPF/Na_Period2/CPF_Period2", title, "_Q4.png", sep=""))
#   polarPlot(data_period2, poll = title, stati = 'cpf', percentile = c(75, 100), main=paste(title, 'Jul17 to Nov14 Q4'))
#   dev.off()
#   print('done Q4')
#   if (title=='Ca'){
#     png(paste("images/CPF/Na_Period2/CPF_Period2", title, "_Q1.png", sep=""))
#     polarPlot(data_period2, poll = title, stati = 'cpf', percentile = c(0, 25), main=paste(title, 'Jul17 to Nov14 Q1'), limits=c(0, 0.025))
#     dev.off()
#   } else if (title=='Co'){
#     png(paste("images/CPF/Na_Period2/CPF_Period2", title, "_Q1.png", sep=""))
#     polarPlot(data_period2, poll = title, stati = 'cpf', percentile = c(0, 25), main=paste(title, 'Jul17 to Nov14 Q1'), limits=c(0, 8.441e-05))
#     dev.off()
#   } else {
#     png(paste("images/CPF/Na_Period2/CPF_Period2", title, "_Q1.png", sep=""))
#     polarPlot(data_period2, poll = title, stati = 'cpf', percentile = c(0, 25), main=paste(title, 'Jul17 to Nov14 Q1'))
#     dev.off()
#   }
#     print('done Q1')
# }
# 
# print(data_period2$Ca[which(data_period2$Ca<=0.0250)])

ggplot(mergeddata)+geom_point(aes(x=date, y=OC_EC, color=Event))

ggplot(mergeddata)+geom_point(aes(x=`C.Orgánico`, y=`C.Elemental`, color=Event))

ggplot(mergeddata)+geom_point(aes(x=`C.Orgánico`, y=nnsK, color=Event))+theme_bw()


ggplot(mergeddata)+geom_line(aes(x=date, y=nssK_EC, color="nssK/EC" ))+
  geom_line(aes(x=date, y=nssK_OC*5, color="nssK/OC"))+
  geom_point(aes(x=date, y=OC_EC/5, color=Event))+
  theme_bw()

ggplot(mergeddata)+geom_point(aes(x=nssK, y=`C.Elemental`, color=Event ))+
  # geom_point(aes(x=date, y=nssK_OC, color="nssK/OC"))+
  # geom_point(aes(x=date, y=OC_EC/10, color="OC/EC"))+
  theme_bw()

ggplot(mergeddata)+geom_point(aes(x=nssK, y=`C.Orgánico`, color="no Event"))+
  geom_point(data=mergeddata[mergeddata$Event %in% c("S","SN","SP","SL","DS"),],
             aes(x=nssK, y=`C.Orgánico`, color="Event"))
  # geom_point(aes(x=date, y=nssK_OC, color="nssK/OC"))+
  # geom_point(aes(x=date, y=OC_EC/10, color="OC/EC"))+
  theme_bw()

ggplot(mergeddata)+geom_point(aes(x=`C.Elemental`, y=`C.Orgánico`, color=Event ))+
  # geom_point(aes(x=date, y=nssK_OC, color="nssK/OC"))+
  # geom_point(aes(x=date, y=OC_EC/10, color="OC/EC"))+
  theme_bw()

ggplot(mergeddata)+
  # geom_point(aes(x=date, y=K, color=Event ))+
  geom_line(aes(x=date, y=K, color="K"))+
  geom_line(aes(x=date, y=nssK, color="nssK"))+
  # geom_point(aes(x=date,y=nssK, color=Event))+
  geom_point(aes(x=date, y=OC_EC/10, color=Event))+
  theme_bw()

ggplot(mergeddata)+geom_point(aes(x=nssK, y=`C.Elemental`, color="no Event" ))+
  geom_point(data=mergeddata[mergeddata$Event %in% c("S","SN","SP","SL"),],
             aes(x=nssK, y=`C.Elemental`, color="Event" ))+
  # geom_point(aes(x=date, y=nssK_OC, color="nssK/OC"))+
  # geom_point(aes(x=date, y=OC_EC/10, color="OC/EC"))+
  theme_bw()

ggplot(mergeddata)+geom_point(aes(x=nssK, y=`C.Orgánico`, color="no Event" ))+
  geom_point(data=mergeddata[mergeddata$Event %in% c("S","SN","SP","SL"),],
             aes(x=nssK, y=`C.Orgánico`, color="Event" ))+
  # geom_point(aes(x=date, y=nssK_OC, color="nssK/OC"))+
  # geom_point(aes(x=date, y=OC_EC/10, color="OC/EC"))+
  theme_bw()

ggplot(mergeddata)+geom_point(aes(x=`C.Elemental`, y=`C.Orgánico`, color="no Event" ))+
  geom_point(data=mergeddata[mergeddata$Event %in% c("S","SN","SP","SL"),],
             aes(x=`C.Elemental`, y=`C.Orgánico`, color="Event" ))+
  # geom_point(aes(x=date, y=nssK_OC, color="nssK/OC"))+
  # geom_point(aes(x=date, y=OC_EC/10, color="OC/EC"))+
  theme_bw()


ggplot(mergeddata)+geom_point(aes(x=NO3, y=NH4, color="no Event" ))+
  geom_point(data=mergeddata[mergeddata$Event %in% c("S","SN","SP","SL"),],
             aes(x=NO3, y=NH4, color="Event" ))+
  # geom_point(aes(x=date, y=nssK_OC, color="nssK/OC"))+
  # geom_point(aes(x=date, y=OC_EC/10, color="OC/EC"))+
  theme_bw()

ggplot(mergeddata)+geom_point(aes(x=SO4, y=NH4, color="no Event" ))+
  geom_point(data=mergeddata[mergeddata$Event %in% c("S","SN","SP","SL"),],
             aes(x=SO4, y=NH4, color="Event" ))+
  # geom_point(aes(x=date, y=nssK_OC, color="nssK/OC"))+
  # geom_point(aes(x=date, y=OC_EC/10, color="OC/EC"))+
  theme_bw()

ggplot(mergeddata)+geom_point(aes(x=Na, y=Cl, color="no Event" ))+
  geom_point(data=mergeddata[mergeddata$Event %in% c("S","SN","SP","SL"),],
             aes(x=Na, y=Cl, color="Event" ))+
  # geom_point(aes(x=date, y=nssK_OC, color="nssK/OC"))+
  # geom_point(aes(x=date, y=OC_EC/10, color="OC/EC"))+
  theme_bw()
# crustal
# "Ca","Al",  "Mg", "Fe", "Ti","Ba", "Mn"
polarPlot(mergeddata, pollutant = "Ca", statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, 
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata, pollutant = "Al", statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed,
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata, pollutant = c("Ca","Al"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, key.footer="[ng/m3]",
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata, pollutant = c( "Fe"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, 
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata, pollutant = c( "Mg"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, 
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata, pollutant = c(  "Ti"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed,
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata, pollutant = c(  "Ba"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, 
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata, pollutant = c(  "Mn"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, 
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata, pollutant = c("V", "Mn","Fe", "Ni", "Cu", "Zn", "Pb", "Ca", "Mg"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, normalise = T, limits = c(0,2),
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata, pollutant = c("V", "Mn","Fe", "Ni", "Cu", "Zn", "Pb", "Ca", "Mg"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, 
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata, pollutant = c("V", "Mn","Fe", "Ni", "Cu", "Zn", "Pb", "Mg"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, 
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata, pollutant = c("K", "NO3","Cl", "SO4", "NH4","Na"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, normalise = T, limits = c(0,2),
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata, pollutant = c("K", "NO3","Cl", "SO4", "NH4","Na"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, 
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata, pollutant = c("K"), type="season", hemisphere="southern",statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed,  k=80, normalise=T,
          mis.col = "transparent", cols=colorset)
ppKJAS=polarPlot(selectByDate(mergeddata, start = "2019-07-01", end = "2019-09-30"), 
          pollutant = c("K"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed,  k=100, 
          mis.col = "transparent", cols=colorset)
ppK=polarPlot(mergeddata, 
                 pollutant = c("K"), statistic = "mean",  min.bin = 2, 
                 upper =upper_windspeed,  k=100, 
                 mis.col = "transparent", cols=colorset)
png(paste0(pathgraphs,"/polarplotppK_",colorset,".png"), width = 717 * 5, height = 339* 5, res = 300)
print(ppK$plot,split=c(1, 1, 2, 1))
print(ppKJAS$plot,split=c(2, 1, 2, 1), newpage=FALSE)
dev.off()

polarPlot(mergeddata, 
              pollutant = c("K","Na"), statistic = "mean",  min.bin = 2, 
              upper =upper_windspeed, normalise = T,
          k=100, limits = c(0,1.5),
              mis.col = "transparent", cols=colorset)
ggplot(selectByDate(mergeddata, start = "2019-07-25", end = "2019-09-25"))+geom_line(aes(x=date, y=K))
corPlot(selectByDate(mergeddata, start = "2019-07-25", end = "2019-09-25"))

polarPlot(mergeddata, pollutant = c("Ca", "Mg"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, normalise = T, limits = c(0,2),
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata, pollutant = c("Cu", "K"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, normalise = T, limits = c(0,2),
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata,pollutant="PM2.5")
polarPlot(selectByDate(mergeddata, start = "2019-07-25", end = "2019-09-25"),
          pollutant="PM2.5")

polarPlot(selectByDate(mergeddata, start = "2019-07-25", end = "2019-09-25"),
          pollutant="NH4")
polarPlot(mergeddata, type="season",pollutant="NH4", hemisphere="southern",min.bin = 2, k=80)
polarPlot(mergeddata, type="season",pollutant="NO3", hemisphere="southern",min.bin = 2, k=80)
polarPlot(mergeddata, type="season",pollutant="SO4", hemisphere="southern",min.bin = 2, k=80)
polarPlot(mergeddata, type="season",pollutant="K", hemisphere="southern",min.bin = 2, k=80)
polarPlot(mergeddata, type="season",pollutant="C.Orgánico", hemisphere="southern",min.bin = 2, k=80)
polarPlot(mergeddata, type="season",pollutant="C.Elemental", hemisphere="southern",min.bin = 2, k=80)
