library(openair)
library(ggplot2)
library(lubridate)
library(readxl)
setwd("~/Documents/paper_facu/Rscripts")
BA_events_testM <- read_excel("../BA_events_testMnew.xlsx")
BA_events_testM$date <- as.POSIXct(BA_events_testM$date, tz='UTC')
BA_events_testM$Event_M <- as.factor(BA_events_testM$Event_M)

factor(BA_events_testM$Event_M)

PMF_BA_full <- read_excel("../data/PMF_BA_fullv3.xlsx",
                          sheet = "CONC", col_types = c("date",
                                                        "numeric", "numeric", "numeric",
                                                        "numeric", "numeric", "numeric",
                                                        "skip", "skip", "numeric", "numeric",
                                                        "numeric", "numeric", "numeric",
                                                        "numeric", "numeric", "numeric",
                                                        "numeric", "numeric", "numeric",
                                                        "numeric", "numeric", "numeric",
                                                        "numeric", "numeric", "numeric",
                                                        "numeric", "numeric", "numeric",
                                                        "numeric", "skip", "skip", "skip",
                                                        "skip", "skip", "skip", "skip",
                                                        "skip", "skip", "skip", "skip", "skip",
                                                        "numeric", "numeric", "numeric"),
                          na ='-999')
# PMF_BA_full0 <- read_excel("PMF_BA_full.xlsx",
#                           sheet = "CONC", col_types = c("date",
#                                                         "numeric", "numeric", "numeric",
#                                                         "numeric", "numeric", "numeric",
#                                                         "skip", "skip", "numeric", "numeric",
#                                                         "numeric", "numeric", "numeric",
#                                                         "numeric", "numeric", "numeric",
#                                                         "numeric", "numeric", "numeric",
#                                                         "numeric", "numeric", "numeric",
#                                                         "numeric", "numeric", "numeric",
#                                                         "numeric", "numeric", "numeric",
#                                                         "numeric", "skip", "skip", "skip",
#                                                         "skip", "skip", "skip", "skip",
#                                                         "skip", "skip", "skip", "skip", "skip",
#                                                         "numeric", "numeric", "numeric"),
#                           na ='-999')
PMF_BA_full=dplyr::rename(PMF_BA_full, EC = 'C Elemental',TC = 'C Total',OC = 'C Orgánico')
PMF_BA_full$KNON=(PMF_BA_full$K-0.6*PMF_BA_full$Fe)
PMF_BA_full$OC_EC = PMF_BA_full$OC/PMF_BA_full$EC
PMF_BA_full$OC_K=PMF_BA_full$OC/PMF_BA_full$K

MWNH4=18.04
MWNO3=62.0049
MWSO4=96.06
MWK=39.0898

PMF_BA_full$NH4molar=PMF_BA_full$NH4/MWNH4
PMF_BA_full$NO3molar=PMF_BA_full$NO3/MWNO3
PMF_BA_full$SO4molar=PMF_BA_full$SO4/MWSO4
PMF_BA_full$Kmolar=PMF_BA_full$K/MWK



PMF_BA_full_so=PMF_BA_full[PMF_BA_full$date<=as.POSIXct('2019-05-23') | PMF_BA_full$date>=as.POSIXct('2019-06-02'),]
pruebacut3_so=selectByDate(PMF_BA_full_so, end="2019-10-01")
PMF_BA_full=merge(PMF_BA_full,BA_events_testM, by="date")

corPlot(PMF_BA_full_so, dendrogram = TRUE,method = "pearson", main= "R pearson, sin outliers")
corPlot(pruebacut[,c("SO4","NH4","NO3","K")], dendrogram = TRUE,method = "pearson", main ="R pearson, matriz completa")
# Iones ####
#Octubre en adelante
pruebacut=selectByDate(PMF_BA_full, start="2019-10-01")
corPlot(pruebacut[,c("SO4","NH4","NO3","K","Na sol","Mg")], dendrogram = TRUE,method = "pearson", main ="desde octubre")
#Hasta Mayo
# 
pruebacut=selectByDate(PMF_BA_full, end="2019-05-20")
corPlot(pruebacut[,c("SO4","NH4","NO3","K","Na sol","Mg")], dendrogram = TRUE, 
        method = "pearson", main ="Hasta Mayo")
ggplot(pruebacut) + geom_point(aes(x=K,y=NO3))

# mayo a octubre
pruebacut=selectByDate(PMF_BA_full, start="2019-05-20",end="2019-10-01")
corPlot(pruebacut[,c("SO4","NH4","NO3","K","Na sol","Mg")], dendrogram = TRUE, 
        method = "pearson", main ="Hasta Mayo")
ggplot(pruebacut) + geom_point(aes(x=`Na sol` ,y=SO4))

pruebacut2=selectByDate(PMF_BA_full, end="2019-10-01")
pruebacut3=selectByDate(PMF_BA_full,  start="2019-05-15", end="2019-10-01")
corPlot(pruebacut3[,c("SO4","NH4","NO3","K")], dendrogram = TRUE,method = "pearson", main ="mayo a octubre")
pruebacut5=selectByDate(PMF_BA_full,  end="2019-05-15")
corPlot(pruebacut5[,c("SO4","NH4","NO3","K")], dendrogram = TRUE,method = "pearson", main ="hasta mayo")

pruebacut6=selectByDate(PMF_BA_full,  start="2019-10-01")
corPlot(pruebacut6[,c("SO4","NH4","NO3","K")], dendrogram = TRUE,method = "pearson", main ="desde octubre")




corPlot(pruebacut2[,c("SO4","NH4","NO3","K")], dendrogram = TRUE,method = "pearson", main ="R pearson, matriz completa")
corPlot(PMF_BA_full_so[,c("SO4","NH4","NO3","K")], dendrogram = TRUE,method = "pearson", main ="R pearson, matriz completa")
corPlot(pruebacut3_so[,c("SO4","NH4","NO3","K")], dendrogram = TRUE,method = "pearson", main ="R pearson, matriz completa")
corPlot(PMF_BA_full[,c("SO4","NH4","NO3","K")], dendrogram = TRUE,method = "pearson", main ="R pearson, matriz completa")
quantile(mydata$`PM2,5`, na.rm = T, 0.9)
quantile(mydata$`PM2,5`, na.rm = T, 0.95)
quantile(mydata$`PM2,5`, na.rm = T, 0.99)
corPlot(PMF_BA_full, dendrogram = TRUE,method = "spearman", main ="R spearman, matriz completa")
corPlot(PMF_BA_full_so, dendrogram = TRUE,method = "spearman", main ="R spearman, sin outliers")

corPlot(PMF_BA_full[,c(2:5,7,10:30,43:46)], dendrogram = TRUE,method = "spearman", main ="R spearman, matriz completa")
corPlot(PMF_BA_full[,c(2:5,7,10:30,43:46)], dendrogram = TRUE,method = "pearson", main ="R pearson, matriz completa")

scatterPlot(PMF_BA_full, x="Ni", y="V")
ggplot(PMF_BA_full)+ geom_point(aes(x=Ni, y=V))+geom_abline(slope = 0.7)+geom_abline(slope = 0.15)+geom_abline(slope =2.8)
# Define event levels
event_levels <- c("SN", "SL", "S", "SC")

# Count the number of events
event_count <- sum(table(BA_events_testM$Event_M[month(BA_events_testM$date)>=3&month(BA_events_testM$date)<=5])[event_levels])
print(event_count)
event_count <- sum(table(BA_events_testM$Event_M[month(BA_events_testM$date)>=6&month(BA_events_testM$date)<=8])[event_levels])
print(event_count)
event_count <- sum(table(BA_events_testM$Event_M[month(BA_events_testM$date)>=9&month(BA_events_testM$date)<=11])[event_levels])
print(event_count)
event_count <- sum(table(BA_events_testM$Event_M[month(BA_events_testM$date)==12|month(BA_events_testM$date)<=2])[event_levels])
print(event_count)

# Traffic related species
corPlot(PMF_BA_full[,c(2,9,16,20,21,26,27,28,29,30,31)], dendrogram = TRUE,method = "pearson", main ="R pearson, matriz completa")
# sacar Mo y Ag


# PMF_BA_full$Lote=as.factor(PMF_BA_full$Lote)
# relaciones ####
ggplot(PMF_BA_full_so[PMF_BA_full_so$Ni!=0,])+geom_point(aes(x=Ni,y=V))+
  geom_abline(aes(slope=2.4, intercept=0), linetype="dashed")+
  geom_abline(aes(slope=3, intercept=0, color="oil combustion from ship engines"), linetype="dashed")+
  geom_abline(aes(slope=1.49, intercept=0, color="road traffic"), linetype="dashed")+
  geom_abline(aes(slope=0.01, intercept=0, color="prueba"), linetype="dashed")

ggplot(PMF_BA_full_so[PMF_BA_full_so$Cd<0.005 & PMF_BA_full_so$Sb>0,])+geom_point(aes(x=Cd,y=Sb))+
  geom_abline(aes(slope=5, intercept=0, color="road brake pad wear "), linetype="dashed")+
  geom_abline(aes(slope=6, intercept=0, color="prueba"), linetype="dashed")

ggplot(PMF_BA_full_so[PMF_BA_full_so$Cd<0.005 & PMF_BA_full_so$Cd>0.00001,],aes(x=Cd,y=Sb))+geom_point()+
  stat_summary(fun.data= mean_cl_normal) + 
  geom_smooth(method='lm',se = F)+geom_abline(aes(slope=0.183260305, intercept=0.002618539))
coefs <- coef(lm(Sb~Cd, data = PMF_BA_full_so[PMF_BA_full_so$Cd>0.00001,]))

# cor plot geological minerals
corPlot(PMF_BA_full_so[,c("Al","Ba", "Ca","Fe", "Mn", "Mg", "Ti","Sb")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson, sin outliers")
corPlot(PMF_BA_full_so[,c( "Ca","Ba",  "Mg", "Ti","Sb")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson, sin outliers")
corPlot(PMF_BA_full_so[,c( "Fe","Al",  "Mn", "Ti")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson, sin outliers")

corPlot(PMF_BA_full_so[,c( "Fe","Al",  "Mn", "Ti")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson, sin outliers")
corPlot(PMF_BA_full_so[,c( "Fe","Al","Ba","Cu","Pb","Sb","Zn", "As")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson, traffic related")
corPlot(PMF_BA_full_so[,c("Zn", "Cu","Pb")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson, traffic related")


# Fe, Ba, Cu, Pb, Sb, Zn, Al
corPlot(PMF_BA_full[,c( "Na sol","Cl",  "NH4", "NO3","SO4", "KNON")], 
        dendrogram = TRUE,method = "spearman", main= "R spearman")
PMF_BA_full$nssK=PMF_BA_full$K-0.6*PMF_BA_full$Fe-0.037*PMF_BA_full$`Na sol`
PMF_BA_full$nssK_OC=PMF_BA_full$nssK/PMF_BA_full$OC
PMF_BA_full$nssK_EC=PMF_BA_full$nssK/PMF_BA_full$EC
PMF_BA_full$OC_EC=PMF_BA_full$OC/PMF_BA_full$EC

PMF_BA_full_so=PMF_BA_full[PMF_BA_full$date<=as.POSIXct('2019-05-23') | PMF_BA_full$date>=as.POSIXct('2019-06-02'),]

corPlot(PMF_BA_full_so[,c( "Na sol","Cl",  "NH4", "NO3","SO4", "nssK")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson")
corPlot(PMF_BA_full_so[,c( "Ca","Al",  "Mg", "Fe", "Ti","Ba", "Mn")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson, sin outliers")


corPlot(PMF_BA_full[,c( "OC","EC","K","NO3","NH4","SO4")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson")
corPlot(PMF_BA_full[PMF_BA_full$Event_F %in% c("SI","SF","SN"),c( "OC","EC","K","NO3","NH4","SO4")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson")
corPlot(PMF_BA_full[PMF_BA_full$Event_F %in% c("no"),c( "OC","EC","K","NO3","NH4","SO4")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson")

corPlot(PMF_BA_full_so[,c( "OC","EC","K","NO3","NH4","SO4")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson")
corPlot(PMF_BA_full_so[PMF_BA_full$Event_F %in% c("SI","SF","SN"),c( "OC","EC","K","NO3","NH4","SO4")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson")
corPlot(PMF_BA_full_so[PMF_BA_full$Event_F %in% c("no"),c( "OC","EC","K","NO3","NH4","SO4")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson")


corPlot(selectByDate(PMF_BA_full, start = "2019-07-25", end = "2019-09-25"),
        pollutants = c("NO3","NH4","K","C Elemental","C Orgánico","SO4"),dendrogram = TRUE)
corPlot(PMF_BA_full,dendrogram = TRUE)

PMF_BA_full$neutralization=PMF_BA_full$NH4/(PMF_BA_full$NO3+PMF_BA_full$SO4)
ggplot(PMF_BA_full)+
  geom_line(aes(x=date, y=neutralization))+
  geom_line(aes(x=date, y=K,color="K"))+
  geom_line(aes(x=date, y=SO4,color="Cu"))
  
meteoobsday=timeAverage(meteoobs,avg.time = "day",statistic = "sum")
TVmeteo=timeVariation(meteoobs,pollutant = "PRECIP 6HS (mm)")
ggplot()+geom_line(data=PMF_BA_full,aes(x=date,y=`PM2,5`))+
  geom_point(data=meteoobsday,aes(x=date,y=`PRECIP 6HS (mm)`*2, color="precip"))
TVmeteoobs=timeVariation(meteoobsday, pollutant = "PRECIP 6HS (mm)", statistic = "median")
print(TVmeteoobs$plot$month)

TVmean=timeVariation(PMF_BA_full,pollutant = "PM2,5")
print(TVmean$plot$month)

timePlot(PMF_BA_full,pollutant = "PM2,5")
timePlot(PMF_BA_full,pollutant = "PM2,5", avg.time = "month")
timePlot(PMF_BA_full,pollutant = "PM2,5", avg.time = "month")

timePlot(PMF_BA_full,pollutant = "OC_K")

# IONES ####
PMF_BA_full$NH4NO3 = (PMF_BA_full$NH4/18.04)/(PMF_BA_full$NO3/62)

png("../images/ratioNH4NO3.png")
ggplot()+geom_line(data=PMF_BA_full,aes(x=date,y=NH4NO3))+geom_abline(slope=0,intercept=1, linetype=3)
dev.off()
PMF_BA_full$NH4SO4 = (PMF_BA_full$NH4/18.04)/(PMF_BA_full$SO4/96.06)
png("../images/ratioNH4SO4.png")
ggplot()+geom_line(data=PMF_BA_full,aes(x=date,y=NH4SO4))+geom_abline(slope=0,intercept=1, linetype=3)
dev.off()

PMF_BA_full$NH4SO4NO3 = (PMF_BA_full$NH4/18.04)/((2*PMF_BA_full$SO4/96.06)+(PMF_BA_full$NO3/62))
png("../images/ratioNH4SO4NO3.png")
ggplot()+geom_line(data=PMF_BA_full,aes(x=date,y=NH4SO4NO3))+geom_abline(slope=0,intercept=1, linetype=3)
dev.off()

PMF_BA_full$NH4SO4NO3 = (PMF_BA_full$NH4/18.04)/((PMF_BA_full$SO4/96.06)+(PMF_BA_full$NO3/62))
png("../images/ratioNH41SO4NO3.png")
ggplot()+geom_line(data=PMF_BA_full,aes(x=date,y=NH4SO4NO3))+geom_abline(slope=0,intercept=1, linetype=3)
dev.off()

PMF_BA_full$KSO4 = (PMF_BA_full$K/39.1)/((PMF_BA_full$SO4/96.06))
ggplot()+geom_line(data=PMF_BA_full,aes(x=date,y=KSO4))+geom_abline(slope=0,intercept=1, linetype=3)
ggplot()+geom_line(data=PMF_BA_full,aes(x=date,y=K))+geom_abline(slope=0,intercept=1, linetype=3)

PMF_BA_full$NO3SO4 = (PMF_BA_full$NO3/62)/((PMF_BA_full$SO4/96.06))
png("../images/ratioNO3SO4.png")
ggplot()+geom_line(data=PMF_BA_full,aes(x=date,y=NO3SO4))+geom_abline(slope=0,intercept=1, linetype=3)
dev.off()

png("../images/ratioSO4.png")
ggplot()+geom_line(data=PMF_BA_full,aes(x=date,y=SO4))
dev.off()
png("images/ratioNO3.png")
ggplot()+geom_line(data=PMF_BA_full,aes(x=date,y=NO3))
dev.off()
png("images/ratioNH4.png")
ggplot()+geom_line(data=PMF_BA_full,aes(x=date,y=NH4))+geom_abline(slope=0,intercept=1, linetype=3)

dev.off()

ggplot(PMF_BA_full)+
  geom_line(aes(x=date, y=Kmolar,color="K"))+
  geom_line(aes(x=date, y=NH4molar,color="NH4"))+
  geom_line(aes(x=date, y=NO3molar,color="NO3"))+
  geom_line(aes(x=date, y=SO4molar,color="SO4"))+
  geom_point(data=PMF_BA_full[PMF_BA_full$Event_F %in% c("SI","SF"),],aes(x=date, y=NH4molar),color="black",show.legend = FALSE)

PMF_BA_full$NH4SO4=PMF_BA_full$NH4molar/PMF_BA_full$SO4molar
PMF_BA_full$NH4NO3=PMF_BA_full$NH4molar/PMF_BA_full$NO3molar
PMF_BA_full$NO3SO4=PMF_BA_full$NO3molar/PMF_BA_full$SO4molar
summary(selectByDate(PMF_BA_full,end = "2019-05-15"))
summary(selectByDate(PMF_BA_full,start = "2019-05-15",end = "2019-10-01"))
summary(selectByDate(PMF_BA_full,start = "2019-10-01"))

ggplot(PMF_BA_full)+
  geom_point(aes(x=NO3molar, y=(NH4molar/SO4molar-2)*SO4molar,color="EXCEDENCIA"))
# hasta aca####
# ggplot()+geom_line(data=PMF_BA_full,aes(x=date,y=K))+geom_abline(slope=0,intercept=1, linetype=3)

PMF_BA_full$bal = (2*(PMF_BA_full$Ca/40.08)+2*(PMF_BA_full$Mg/24.31)+(PMF_BA_full$`Na sol`/23)+(PMF_BA_full$K/39.1)+(PMF_BA_full$NH4/18.04))/((2*PMF_BA_full$SO4/96.06)+(PMF_BA_full$NO3/62)+(PMF_BA_full$Cl/35.45))
ggplot()+geom_line(data=PMF_BA_full,aes(x=date,y=bal))+geom_abline(slope=0,intercept=1, linetype=3)

colMeans(PMF_BA_full[,2:40],na.rm = T)

PMF_BA_full$bal = ((PMF_BA_full$`Na sol`/23)+(PMF_BA_full$K/39.1)+(PMF_BA_full$NH4/18.04))/((2*PMF_BA_full$SO4/96.06)+(PMF_BA_full$NO3/62)+(PMF_BA_full$Cl/35.45))
ggplot()+geom_line(data=PMF_BA_full,aes(x=date,y=bal))+geom_abline(slope=0,intercept=1, linetype=3)

ggplot()+geom_line(data=PMF_BA_full,aes(x=date,y=SO4))+geom_abline(slope=0,intercept=1, linetype=3)

# 
gases20192020day=timeAverage(gases20192020, avg.time = "day")

gases20192020$date=gases20192020$day+gases20192020$Hora
ggplot(selectByDate(PMF_BA_full,start = "2019-05-15",end = "2019-06-15"))+
  geom_line(aes(x=date,y=SO4,color="SO4"))+
  geom_line(aes(x=date,y=NO3,color="NO3"))+
  geom_line(aes(x=date,y=`Na sol`,color="Na"))+
  geom_line(aes(x=date,y=Cl,color="Cl"))+
  geom_line(aes(x=date,y=NH4,color="NH4"))

ggplot(PMF_BA_full)+geom_line(aes(x=date,y=SO4,color="SO4"))+geom_point(aes(x=date,y=SO4,color=OrigenTag))
ggplot(PMF_BA_full)+geom_line(aes(x=date,y=NO3,color="NO3"))+geom_point(aes(x=date,y=NO3,color=OrigenTag))
ggplot(PMF_BA_full)+geom_line(aes(x=date,y=EC,color="EC"))+geom_point(aes(x=date,y=EC,color=OrigenTag))
