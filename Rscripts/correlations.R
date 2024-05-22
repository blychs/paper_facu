library(openair)
library(ggplot2)
library(lubridate)
library(readxl)
BA_events_testM <- read_excel("BA_events_testM.xlsx")
BA_events_testM$date <- as.POSIXct(BA_events_testM$date, tz='UTC')
BA_events_testM$Event_M <- as.factor(BA_events_testM$Event_M)

factor(BA_events_testM$Event_M)

PMF_BA_full <- read_excel("data/PMF_BA_fullv3.xlsx",
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
PMF_BA_full0 <- read_excel("PMF_BA_full.xlsx",
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
PMF_BA_full=dplyr::rename(PMF_BA_full, EC = 'C Elemental',TC = 'C Total',OC = 'C Orgánico')
PMF_BA_full$KNON=(PMF_BA_full$K-0.6*PMF_BA_full$Fe)
PMF_BA_full$OC_EC = PMF_BA_full$OC/PMF_BA_full$EC
PMF_BA_full$OC_K=PMF_BA_full$OC/PMF_BA_full$K

PMF_BA_full_so=PMF_BA_full[PMF_BA_full$date<=as.POSIXct('2019-05-23') | PMF_BA_full$date>=as.POSIXct('2019-06-02'),]
corPlot(PMF_BA_full_so, dendrogram = TRUE,method = "pearson", main= "R pearson, sin outliers")
corPlot(PMF_BA_full, dendrogram = TRUE,method = "pearson", main ="R pearson, matriz completa")

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

PMF_BA_full=merge(PMF_BA_full,BA_events_testM, by="date")
PMF_BA_full$Lote=as.factor(PMF_BA_full$Lote)
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
