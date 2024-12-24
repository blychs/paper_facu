library(openair)
# 00 Carga de datos ####
library(ggplot2)
library(lubridate)
library(readxl)
setwd("~/mdiaz/Documents/paper_facu")
pathgraphs="Figures"
upper_windspeed=8.5
colorset="plasma"
# En data_every_hour_obs.csv solo estan actualizados los datos de meteo el resto del archivo es viejo
data <- read.csv("data_every_hour_obsv5.csv")                 
data$date <- as.POSIXct(data$date, tz='UTC')
data$temp[data$temp>900]=NA
data$day<- format(data$date-12*3600, "%Y-%m-%d UTC")
data$day <- as.POSIXct(data$day, tz='UTC')
mergeddata <- data[,-c(8,9,31:42)]
mergeddata$Sb_ng=mergeddata$Sb*1000
mergeddata$As_ng=mergeddata$As*1000
mergeddata$nssK=mergeddata$K-0.6*mergeddata$Fe-0.037*mergeddata$Na
mergeddata$nssK_OC=mergeddata$nssK/mergeddata$C.Orgánico
mergeddata$nssK_EC=mergeddata$nssK/mergeddata$C.Elemental
mergeddata$OC_EC=mergeddata$C.Orgánico/mergeddata$C.Elemental
mergeddata=dplyr::rename(mergeddata,Event= Event_F)
mergeddata=dplyr::rename(mergeddata,Na= Na.sol)
BA_events_testM <- read_excel("BA_events_testMnew2.xlsx")
BA_events_testM$date <- as.POSIXct(BA_events_testM$date, tz='UTC')


datamerged=merge(mergeddata,BA_events_testM,by.x="day", by.y ="date")
mergeddata$OC_EC = mergeddata$`C.Orgánico`/mergeddata$`C.Elemental`
mergeddata$K_OC=mergeddata$K/mergeddata$`C.Orgánico`
# datamergedcut=datamerged[,c(2:29,42:44,46:51,59:63)]

mergeddata$lat= -34.5730#-0.01
mergeddata$lon=-58.5127#-0.01
mergeddata$Vng=mergeddata$V*1000
# polar maps####
library(openairmaps)
library(leaflet)
customIcon <- makeIcon(
  iconUrl = "/home/usuario/mdiaz/Documents/paper_facu/central.png",  # Cambia esto a la ruta de tu PNG
  iconWidth = 120,  # Ajusta el tamaño según sea necesario
  iconHeight = 120
)

customIconships <- makeIcon(
  iconUrl = "/home/usuario/mdiaz/Documents/paper_facu/ship-flat-boat-by-Vexels.svg",  # Cambia esto a la ruta de tu PNG
  iconWidth = 120,  # Ajusta el tamaño según sea necesario
  iconHeight = 120
)

customIcontree <- makeIcon(
  iconUrl = "/home/usuario/mdiaz/Documents/paper_facu/tree.png",  # Cambia esto a la ruta de tu PNG
  iconWidth = 120,  # Ajusta el tamaño según sea necesario
  iconHeight = 120
)


customIconcars <- makeIcon(
  iconUrl = "/home/usuario/mdiaz/Documents/paper_facu/car-fleet-12792.png",  # Cambia esto a la ruta de tu PNG
  iconWidth = 120,  # Ajusta el tamaño según sea necesario
  iconHeight = 120
)

leaflet(data = mergeddata) %>%
  addTiles() %>%
  addProviderTiles(providers$OpenStreetMap) %>% 
  addPolarMarkers("Vng", 
                  fun = openair::polarPlot,
                  group = "Polar Plot",
                  cols="inferno",
                  alpha = 1,
                  key = FALSE,
                  key.position="left",
                  key.footer="",
                  key.header = "V [ng/m3]"
  )%>%
  addMarkers(lng = -58.344426375249924, lat = -34.64608663021544, 
             popup = "Central Costanera", icon = customIcon) %>%
  addMarkers(lng = -58.380487846322495, lat = -34.57504109311285, 
             popup = "Central Puerto", icon = customIcon )%>%
  addMarkers(lng = -58.37137834376634, lat = -34.51275959413225, 
             popup = "Ships", icon = customIconships ) %>%
  addMarkers(lng = -58.39237834376634, lat = -34.51075959413225, 
             popup = "Ships", icon = customIconships ) %>%  
  addMarkers(lng = -58.38137834376634, lat = -34.518275959413225, 
            popup = "Ships", icon = customIconships ) %>%
  addMarkers(lng = -58.48137834376634, lat = -34.54075959413225, 
             popup = "Car", icon = customIconcars )

polarPlot(mergeddata,pollutant = "Vng",   cols="inferno",        key.position="left",
          key.footer="",
          key.header = "V [ng/m3]")

# Generar el mapa de trayectorias ucustomIconcars# Generar el mapa de trayectorias usando trajMap
trajMap(selectByDate(subset(traj500, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat), 
                     start = "2019-08-26", end = "2019-08-30"),
        origin = TRUE,  
        grid.col = "transparent", 
        map.cols = "transparent",
        projection = "stereographic", 
        orientation = c(0, -65, 0), 
        parameters = NULL)

# Superponer el ícono personalizado con leaflet
leafletProxy("map") %>%
  addMarkers(lng = -58.344426375249924, lat = -34.64608663021544, 
             popup = "Central Costanera", icon = customIcon) %>%
  addMarkers(lng = -58.380487846322495, lat = -34.57504109311285, 
             popup = "Central Puerto", icon = customIcon)

trajMap(selectByDate(subset(traj500, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat), 
                       start = "2019-08-26", end = "2019-08-30"), origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                  projection = "stereographic",   orientation=c(0,-65,0), parameters = NULL)

polarMap(
  mergeddata,
  pollutant = "V",
  x = "ws",
  limits = "free",
  upper = "fixed",
  crs = 4326,
  # type = NULL,
  # popup = NULL,
  # label = NULL,
  provider = "OpenStreetMap",
  cols="inferno",
  alpha = 1,
  key = TRUE,
  key.footer="[ug/m3]",
  key.header = "V",
  legend = TRUE,
  # legend.position = NULL,
  # legend.title = NULL,
  legend.title.autotext = TRUE,
  control.collapsed = TRUE,
  control.position = "topright",
  control.autotext = TRUE,
  d.icon = 200,
  d.fig = 3.5,
  static = FALSE,
  static.nrow = NULL,
  progress = TRUE
)

# Dock sud -34.648961551867856, -58.34236555984067
# Central costanera -34.64608663021544, -58.344426375249924
# Central Puerto -34.57504109311285, -58.380487846322495

# 01 EC OC ####
PPPM25<-polarPlot(mergeddata, pollutant = "PM2.5", statistic = "mean",  min.bin = 2, 
                  upper =upper_windspeed, key.footer="[ug/m3]",key.header = "PM2.5",mis.col = "transparent",
                  cols="inferno")
polarPlot(mergeddata, pollutant = "OC_EC", statistic = "mean",  min.bin = 2, main='OC/EC')
ppEC<-polarPlot(mergeddata, pollutant = "C.Elemental", statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, key.footer="[ug/m3]",key.header = "EC",mis.col = "transparent",
          cols=colorset)
ppOC<-polarPlot(mergeddata, pollutant = "C.Orgánico", statistic = "mean",  min.bin = 2, 
                upper =upper_windspeed, key.footer="[ug/m3]",key.header = "OC",mis.col = "transparent",
                cols=colorset)
# ppOC<-polarPlot(mergeddata, pollutant = "OC_EC", statistic = "mean",  min.bin = 2, 
#                 upper =upper_windspeed, key.footer="[ug/m3]",key.header = "OC/EC",mis.col = "transparent",
#                 cols=colorset)
polarPlot(mergeddata, pollutant = "K_OC", statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, key.footer="",key.header = "K/OC",mis.col = "transparent",
          cols=colorset)
# ppOM<-polarPlot(mergeddata, pollutant = "OM", statistic = "mean",  min.bin = 2, 
#                 upper =upper_windspeed, key.footer="[ug/m3]",key.header = "OM",mis.col = "transparent",
#                 cols=colorset)
png(paste0(pathgraphs,"/polarplotECOCenmug_",colorset,".png"), width = 717 * 5, height = 339* 5, res = 300)
print(ppEC$plot,split=c(1, 1, 2, 1))
print(ppOC$plot,split=c(2, 1, 2, 1), newpage=FALSE)
dev.off()
png(paste0(pathgraphs,"/polarplotSO4_",colorset,".png"), width =350 * 4, height = 350* 4, res = 300)
polarPlot(mergeddata, pollutant = "SO4", statistic = "mean",  min.bin = 2, 
                  upper =upper_windspeed, key.footer="[ug/m3]",key.header = "SO4",mis.col = "transparent",
                  cols=colorset,auto.text = TRUE)
dev.off()
png(paste0(pathgraphs,"/polarplotV_",colorset,".png"), width =350 * 4, height = 350* 4, res = 300)
polarPlot(mergeddata, pollutant = "V", statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, key.footer="[ug/m3]",key.header = "V",mis.col = "transparent",
          cols=colorset,auto.text = TRUE)
dev.off()
# png(paste0(pathgraphs,"/polarplotECOCOM_",colorset,".png"), width = 900 * 4, height = 270* 4, res = 300)
# print(ppEC$plot,split=c(1, 1, 3, 1))
# print(ppOC$plot,split=c(2, 1, 3, 1), newpage=FALSE)
# print(ppOM$plot,split=c(3, 1, 3, 1), newpage=FALSE)
# dev.off()
# 02 As y Sb ####
ppAs<-polarPlot(mergeddata, pollutant = "As_ng", statistic = "mean",  min.bin = 2, 
                upper =upper_windspeed, key.header="As",key.footer="[ng/m3]",
                mis.col = "transparent", cols=colorset)
ppSb<-polarPlot(mergeddata, pollutant = "Sb_ng", statistic = "mean",  min.bin = 2, 
                upper =upper_windspeed, key.header ="Sb",key.footer="[ng/m3]",
                mis.col = "transparent", cols=colorset)
png(paste0(pathgraphs,"/polarplotAsSbng_",colorset,".png"), width = 717 * 4, height = 339* 4, res = 300)
print(ppAs$plot,split=c(1, 1, 2, 1))
print(ppSb$plot,split=c(2, 1, 2, 1), newpage=FALSE)
dev.off()

# 02 Zn, Cu y Pb ####
ppZn<-polarPlot(mergeddata, pollutant = "Zn", statistic = "mean",  min.bin = 2, 
                upper =upper_windspeed, key.header="Zn",key.footer="[ug/m3]",
                mis.col = "transparent", cols=colorset)
ppPb<-polarPlot(mergeddata, pollutant = "Pb", statistic = "mean",  min.bin = 2, 
                upper =upper_windspeed, key.header="Pb",key.footer="[ug/m3]",
                mis.col = "transparent", cols=colorset)
ppCu<-polarPlot(mergeddata, pollutant = "Cu", statistic = "mean",  min.bin = 2, 
                upper =upper_windspeed, key.header="Cu",key.footer="[ug/m3]",
                mis.col = "transparent", cols=colorset)
png(paste0(pathgraphs,"/polarplotPbZnCu_",colorset,".png"), width = 717 * 5, height = 339* 5, res = 300)
print(ppZn$plot,split=c(1, 1, 3, 1))
print(ppPb$plot,split=c(2, 1, 3, 1), newpage=FALSE)
print(ppCu$plot,split=c(3, 1, 3, 1), newpage=FALSE)
dev.off()

# 03 Ni V ####
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


# 04 K ####
ppK<-polarPlot(mergeddata, pollutant = "K", statistic = "mean",  min.bin = 2, 
               upper =upper_windspeed, key.header="K",key.footer="[ug/m3]",
               mis.col = "transparent", cols=colorset)

png(paste0(pathgraphs,"/polarplotK_",colorset,".png"), width = 500 * 5, height = 500* 5, res = 300)
print(ppK$plot)
dev.off()

# 05 graficos auxiliares con datos horarios ####
ggplot(mergeddata)+geom_point(aes(x=date, y=OC_EC, color=Event))

ggplot(mergeddata)+geom_point(aes(x=`C.Orgánico`, y=`C.Elemental`, color=Event))


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
  geom_point(data=mergeddata[mergeddata$Event %in% c("SI","SF"),],
             aes(x=NO3, y=NH4, color="Event" ))+
  # geom_point(aes(x=date, y=nssK_OC, color="nssK/OC"))+
  # geom_point(aes(x=date, y=OC_EC/10, color="OC/EC"))+
  theme_bw()

ggplot(mergeddata)+geom_point(aes(x=SO4, y=NH4, color="no Event" ))+
  geom_point(data=mergeddata[mergeddata$Event %in% c("SI","SF"),],
             aes(x=SO4, y=NH4, color="Event" ))+
  # geom_point(aes(x=date, y=nssK_OC, color="nssK/OC"))+
  # geom_point(aes(x=date, y=OC_EC/10, color="OC/EC"))+
  theme_bw()

ggplot(mergeddata)+geom_point(aes(x=Na, y=Cl, color="no Event" ))+
  geom_point(data=mergeddata[mergeddata$Event %in% c("SI","SF"),],
             aes(x=Na, y=Cl, color="Event" ))+
  theme_bw()
# 06 Polar plots crustal ####
# "Ca","Al",  "Mg", "Fe", "Ti","Ba", "Mn"

ggplot(PMF_BA_full)+geom_line(aes(x=date,y=Al,color="Al"))+
  geom_line(aes(x=date,y=Fe,color="Fe"))+
  geom_line(aes(x=date,y=Ca,color="Ca"))+
  geom_line(aes(x=date,y=Mg,color="Mg"))
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
polarPlot(mergeddata, pollutant = c(  "K"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed,
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata, pollutant = c(  "Ba"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, 
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata, pollutant = c(  "Sb"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, 
          mis.col = "transparent", cols=colorset)
png(paste0(pathgraphs,"/polarplotvariosmetalesnorm_",colorset,".png"), width = 717 * 5, height = 339* 5, res = 300)
polarPlot(mergeddata, pollutant = c("V", "Mn","Fe", "Ni", "Cu", "Zn", "Pb", "Ca", "Mg"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, normalise = T, limits = c(0,2),
          mis.col = "transparent", cols=colorset)
dev.off()

polarPlot(mergeddata, pollutant = c("V", "Mn","Fe", "Ni", "Cu", "Zn", "Pb", "Ca", "Mg"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, 
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata, pollutant = c("V", "Mn","Fe", "Ni", "Cu", "Zn", "Pb", "Mg"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, 
          mis.col = "transparent", cols=colorset)

polarPlot(mergeddata, pollutant = c( "Cu", "Zn", "Pb", "Sb"), 
          statistic = "mean",  min.bin = 2, normalise = T,
          upper =upper_windspeed, mis.col = "transparent", 
          cols=colorset)
polarPlot(mergeddata, pollutant = c( "Fe", "Ba", "Pb"), 
          statistic = "mean",  min.bin = 2, normalise = T,
          upper =upper_windspeed, mis.col = "transparent", 
          cols=colorset)
ggplot(PMF_BA_full)+geom_line(aes(x=date,y=Cu/2,color="Cu"))+
  geom_line(aes(x=date,y=Zn/5,color="Zn"))+geom_line(aes(x=date,y=K/100,color="K"))+
  geom_line(aes(x=date,y=3*Sb,color="Sb"))


ppiones <- polarPlot(mergeddata, pollutant = c("K", "NO3","Cl", "SO4", "NH4","Na"), statistic = "mean",  min.bin = 2, 
           upper =upper_windspeed, normalise = T, limits = c(0,2),
           mis.col = "transparent", cols=colorset)

png(paste0(pathgraphs,"/polarplotions_",colorset,".png"), width = 500 * 5, height = 500* 5, res = 300)
print(ppiones$plot)
dev.off()

polarPlot(mergeddata, pollutant = c("K", "NO3","Cl", "SO4", "NH4","Na"), statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed, 
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata, pollutant = c("K"), type="season", hemisphere="southern",statistic = "mean",  min.bin = 2, 
          upper =upper_windspeed,  k=80, normalise=T,
          mis.col = "transparent", cols=colorset)
polarPlot(mergeddata, pollutant = c("nssK_OC"), type="season", hemisphere="southern",statistic = "mean",  min.bin = 2, 
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
polarPlot(mergeddata, type="season",pollutant="Cl", hemisphere="southern",min.bin = 2, k=80)
polarPlot(mergeddata, type="season",pollutant="K", hemisphere="southern",min.bin = 2, k=80)
polarPlot(mergeddata, type="season",pollutant="C.Orgánico", hemisphere="southern",min.bin = 2, k=80)
polarPlot(mergeddata, type="season",pollutant="C.Elemental", hemisphere="southern",min.bin = 2, k=80)

polarPlot(selectByDate(mergeddata, start = "2019-05-20", end = "2019-06-30"),
          pollutant="NH4")
polarPlot(selectByDate(mergeddata, start = "2019-05-20", end = "2019-06-30"),
          pollutant="NO3")
polarPlot(selectByDate(mergeddata, start = "2019-05-20", end = "2019-06-30"),
          pollutant="SO4")
polarPlot(selectByDate(mergeddata, start = "2019-05-20", end = "2019-06-30"),
          pollutant="Cl")

polarPlot(mergeddata,pollutant="NH4", type="season",hemisphere="southern",min.bin = 2, k=80)
polarPlot(mergeddata, pollutant="NO3",type="season", hemisphere="southern",min.bin = 2, k=80)
polarPlot(mergeddata, pollutant="SO4", type="season",hemisphere="southern",min.bin = 2, k=80)
