library(readr)
library(openair)
library(ggplot2)
library(lubridate)
library(readxl)
setwd("~/mdiaz/Documents/paper_facu/Rscripts")
Simon_mass <- read_csv("../Simon_2011_events_mass.csv",
                       col_types = cols(trace_elements = col_skip()))

Simon_mass$date <- as.POSIXct(Simon_mass$date, tz='UTC')

# Lichtig_mass <- read_csv("../Lichtig_2024_events_mass.csv",
#                        col_types = cols(trace_elements = col_skip()))
# 
# Lichtig_mass$date <- as.POSIXct(Lichtig_mass$date, tz='UTC')


BA_events_testM <- read_excel("../BA_events_testMnew2.xlsx")
BA_events_testM$date <- as.POSIXct(BA_events_testM$date, tz='UTC')
BA_events_testM$Event_F <- as.factor(BA_events_testM$Event_F)

factor(BA_events_testM$Event_F)

PMF_BA_full <- read_excel("../data/PMF_BA_fullv4.xlsx",
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

Simon_mass$PM2.5=rowSums(Simon_mass[,2:8])
summary(Simon_mass)
# Lichtig_mass$PM2.5=rowSums(Lichtig_mass[,2:8])
# Lichtig_mass$diff=Simon_mass$PM2.5-Lichtig_mass$PM2.5
ggplot(Lichtig_mass)+geom_line(aes(x=date,y=inorganic_ions))+
  geom_point(data=BA_events_testM,aes(x=date, y=Simon_mass$inorganic_ions, color =OrigenTag))+
  geom_line(data=Simon_mass[BA_events_testM$Event_F %in% c("SI","SF")], aes(x=date, y=inorganic_ions, color="Simon"))
 
ggplot(Lichtig_mass)+geom_line(aes(x=date, y=diff))
df=merge(Simon_mass,PMF_BA_full,by="date")
df=merge(df,BA_events_testM,by="date")
modStats(df, mod="PM2.5", obs = "PM2,5" )

summary(df)
df$diff=df$PM2.5-df$`PM2,5`
df$diff_perc=df$diff/df$`PM2,5`*100
df$NH4SO4=(df$NH4/MWNH4)/(df$SO4/MWSO4)
ggplot(df)+geom_line(aes(x=date,y=diff))+geom_point(aes(x=date,y=diff,color=NH4SO4))

ggplot(df)+geom_line(aes(x=date,y=diff))
ggplot(df)+geom_line(aes(x=date,y=diff/`PM2,5`*100))

sd(df$PM2.5,na.rm = T)

df$ii_perc=df$inorganic_ions/df$PM2.5*100
df$om_perc=df$organic_mass/df$PM2.5*100
df$gm_perc=df$geological_minerals/df$PM2.5*100
df$ec_perc=df$elemental_C/df$PM2.5*100
df$ss_perc=df$salt/df$PM2.5*100
df$r_perc=df$residual/df$PM2.5*100
# df=df[complete.cases(df),]
df$OCEC=df$`C Orgánico`/df$`C Elemental`
summary(df)
summary(df[df$Event_F %in% c("SI","SF","SO"),]) # OM 13.88 # 5.6 OC
summary(df[df$Event_F %in% c("no"),]) # OM 8.91 # 4.86 OC
summary(selectByDate(df, start = "2019-12-01",end="2020-04-01")) #0.9177
summary(selectByDate(df, start = "2019-06-01",end="2019-09-01")) #0.92
summary(selectByDate(df, start = "2019-04-01",end="2019-06-01")) #1.2
summary(selectByDate(df, start = "2019-09-01",end="2019-12-01")) #0.92



timeVariation(df, pollutant = "organic_mass")
Simon_2011_original_mass <- read_csv("Documents/paper_facu/Simon_2011_original_mass.csv")
Simon_2011_original_mass$date <- as.POSIXct(Simon_2011_original_mass$date, tz='UTC')
Simon_2011_original_mass$PM2.5=rowSums(Simon_2011_original_mass[,2:9])
df=merge(Simon_2011_original_mass,PMF_BA_full,by="date")
modStats(df, mod="PM2.5", obs = "PM2,5" )


Simon_2011_original_mass <- read_csv("Documents/paper_facu/Simon_2011_alltogether_mass.csv")
Simon_2011_original_mass$date <- as.POSIXct(Simon_2011_original_mass$date, tz='UTC')
Simon_2011_original_mass$PM2.5=rowSums(Simon_2011_original_mass[,2:9])
df=merge(Simon_2011_original_mass,PMF_BA_full,by="date")
modStats(df, mod="PM2.5", obs = "PM2,5" )

Simon_2011_original_mass <- read_csv("Documents/paper_facu/Simon_2011_events_mass.csv")
Simon_2011_original_mass$date <- as.POSIXct(Simon_2011_original_mass$date, tz='UTC')
Simon_2011_original_mass$PM2.5=rowSums(Simon_2011_original_mass[,2:9])
df=merge(Simon_2011_original_mass,PMF_BA_full,by="date")
modStats(df, mod="PM2.5", obs = "PM2,5" )
ggplot(df,aes(x=PM2.5, y=`PM2,5`))+geom_point()+geom_abline(slope = 1,intercept = 0,linetype=3)+
xlim(0,50)+ylim(0,50)

Simon_2011_original_mass <- read_csv("Documents/paper_facu/Malm_1994_original_mass.csv")
Simon_2011_original_mass$date <- as.POSIXct(Simon_2011_original_mass$date, tz='UTC')
Simon_2011_original_mass$PM2.5=rowSums(Simon_2011_original_mass[,2:9])
ggplot(df,aes(x=PM2.5, y=`PM2,5`))+geom_point()+geom_abline(slope = 1,intercept = 0,linetype=3)+
xlim(0,50)+ylim(0,50)
df=merge(Simon_2011_original_mass,PMF_BA_full,by="date")
modStats(df, mod="PM2.5", obs = "PM2,5" )


Simon_2011_original_mass <- read_csv("Documents/paper_facu/Malm_1994_original_mass.csv")
Simon_2011_original_mass$date <- as.POSIXct(Simon_2011_original_mass$date, tz='UTC')
Simon_2011_original_mass$PM2.5=rowSums(Simon_2011_original_mass[,2:9])
ggplot(df,aes(x=PM2.5, y=`PM2,5`))+geom_point()+geom_abline(slope = 1,intercept = 0,linetype=3)+
  xlim(0,50)+ylim(0,50)
df=merge(Simon_2011_original_mass,PMF_BA_full,by="date")
modStats(df, mod="PM2.5", obs = "PM2,5" )
ggplot(df)+geom_line(aes(x=date,y=PM2.5, color="Simon_2011"))+geom_line(aes(x=date,y=`PM2,5`, color="PM2.5"))+geom_point(aes(x=date,y=`PM2,5`, color=Origen))
ggplot(df)+geom_line(aes(x=date,y=PM2.5, color="Simon_2011"))+geom_line(aes(x=date,y=`PM2,5`, color="PM2.5"))+geom_point(aes(x=date,y=`PM2,5`, color=Event_F))
ggplot(df)+geom_line(aes(x=date,y=organic_mass, color="PM2.5"))+geom_point(aes(x=date,y=organic_mass, color=Tag))


ggplot(df)+geom_line(aes(x=date,y=NO3, color="NO3"))+geom_line(aes(x=date,y=SO4, color="SO4"))+geom_line(aes(x=date,y=NH4, color="NH4"))+
  geom_line(aes(x=date,y=Cl, color="Cl"))+geom_line(aes(x=date,y=`Na sol`, color="Na"))+geom_line(aes(x=date,y=`C Orgánico`/5, color="OC"))

ggplot(PMF_BA_full)+geom_line(aes(x=date,y=NO3, color="NO3"))+geom_line(aes(x=date,y=SO4, color="SO4"))+geom_line(aes(x=date,y=NH4, color="NH4"))+
  geom_line(aes(x=date,y=Cl, color="Cl"))+geom_line(aes(x=date,y=`Na sol`, color="Na"))+
  geom_point(data=PMF_BA_full[PMF_BA_full$Event_F %in% c("SI","SF"),],aes(x=date,y=NO3, color=Event_F))


ggplot(PMF_BA_full)+geom_line(aes(x=date,y=NO3, color="NO3"))+geom_line(aes(x=date,y=SO4, color="SO4"))+
  geom_point(data=PMF_BA_full[PMF_BA_full$Event_F %in% c("SI","SF"),],aes(x=date,y=NO3, color=OrigenTag))

ggplot(PMF_BA_full)+geom_line(aes(x=date,y=SO4+NH4+NO3,color="Lichtig_mass"))+ geom_line(aes(x=date,y=1.375*SO4+1.29*NO3,color="Simon_mass"))+ 
  geom_line(data=Simon_mass,aes(x=date,y=inorganic_ions,color="Simon"))+
  geom_line(data=Lichtig_mass,aes(x=date,y=inorganic_ions,color="Lichtig"))
  
