# load libraries ####

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
# library(tidyverse)
library(patchwork)
library(readr)
library(openair)
library(lubridate)
library(readxl)
library(openairmaps)
setwd("~/mdiaz/Documents/paper_facu/Rscripts")
source('opentrajutils.R')

# load dataframes ####
Simon_mass <- read_csv("../methods_csvs/Simon_2011_events_mass.csv",
                       col_types = cols(trace_elements = col_skip()))
Simon_mass$date <- as.POSIXct(Simon_mass$date, tz='UTC')
Simon_mass$PM2.5=rowSums(Simon_mass[,2:8])
Simon_mass$ii_perc=Simon_mass$inorganic_ions/Simon_mass$PM2.5*100
Simon_mass$om_perc=Simon_mass$organic_mass/Simon_mass$PM2.5*100
Simon_mass$gm_perc=Simon_mass$geological_minerals/Simon_mass$PM2.5*100
Simon_mass$ec_perc=Simon_mass$elemental_C/Simon_mass$PM2.5*100
Simon_mass$ss_perc=Simon_mass$salt/Simon_mass$PM2.5*100
Simon_mass$others_perc=Simon_mass$others/Simon_mass$PM2.5*100
Simon_mass$r_perc=Simon_mass$residual/Simon_mass$PM2.5*100

BA_events_testM <- read_excel("../data/BA_events_testMnew2.xlsx")
BA_events_testM$date <- as.POSIXct(BA_events_testM$date, tz='UTC')
BA_events_testM$Event_F <- as.factor(BA_events_testM$Event_F)

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
                                                        "numeric", "numeric", "numeric"), na ='-999')
PMF_BA_full$date <- as.POSIXct(PMF_BA_full$date, tz='UTC')
PMF_BA_full$color <- ifelse(BA_events_testM$Event_F %in% c("SI", "SF","SO"), "BB samples", "non-BB samples")
PMF_BA_full$datefactor=factor(PMF_BA_full$date)
PMF_BA_full$date2 <- as.Date(PMF_BA_full$date, format = "%Y-%m-%d")

pathgraphs="Figures"
upper_windspeed=8.5
pscflims = c(0,1)
colorset="plasma"
pathtraj="../tdumps"
pattern ="tdump*"

traj500<-read.tdumps(pathtraj = pathtraj, pattern=pattern)
traj500=traj500[traj500$receptor==3,]

fires <- read_csv("/media/usuario/32b62ac8-a81b-4630-bf83-822c71ee0cad/mdiaz/Documents/paper_facu/data/fires/DL_FIRE_M-C61_563912/fire_archive_M-C61_563912.csv")
fires$acq_time <- sprintf("%04d", as.integer(fires$acq_time))  
fires$hour <- substr(fires$acq_time, 1, 2)  
fires$minute <- substr(fires$acq_time, 3, 4)  
fires$datetime <- as.POSIXct(paste(fires$acq_date, fires$hour, fires$minute), 
                             format = "%Y-%m-%d %H %M")

start_date <- as.POSIXct("2019-08-26")
end_date <- as.POSIXct("2019-08-30")
fires$datetime <- as.POSIXct(fires$datetime)

# Crea el rango temporal de fuegos a plotear
start_filter <- start_date - as.difftime(72, units = "hours")  # 72 hs antes
end_filter <- end_date   # 24 hs después
fires_filtered <- fires %>%
  filter(datetime >= start_filter & datetime <= end_filter)

pathgraphs="Figures"
upper_windspeed=8.5
colorset="plasma"
# En data_every_hour_obs.csv solo estan actualizados los datos de meteo el resto del archivo es viejo
data <- read.csv("../data/data_every_hour_obsv5.csv")                 
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
datamerged=merge(mergeddata,BA_events_testM,by.x="day", by.y ="date")
mergeddata$OC_EC = mergeddata$`C.Orgánico`/mergeddata$`C.Elemental`
mergeddata$K_OC=mergeddata$K/mergeddata$`C.Orgánico`
# datamergedcut=datamerged[,c(2:29,42:44,46:51,59:63)]

mergeddata$lat= -34.5730#-0.01
mergeddata$lon=-58.5127#-0.01
mergeddata$Vng=mergeddata$V*1000


# prepare data ####
# Simon_mass$residual[Simon_mass$residual<0]=0
# mean_data <- Simon_mass %>%
#   mutate(Event_F = PMF_BA_full$color) %>%
#   summarise(across(c(ii_perc, om_perc, ec_perc, gm_perc, ss_perc, others_perc, r_perc),\(x) mean(x, na.rm = TRUE))) %>%
#   pivot_longer(cols = c(ii_perc, om_perc, ec_perc, gm_perc, ss_perc, others_perc, r_perc),
#                names_to = "category", values_to = "value") %>%
#   mutate(category = recode(category,
#                            "ii_perc" = "Inorganic Ions \n 10%",
#                            "om_perc" = "Organic Matter \n 57%",
#                            "ec_perc" = "Elemental Carbon \n 6%",
#                            "gm_perc" = "Geological Minerals \n 14%",
#                            "ss_perc" = "Sea Salt \n 3%",
#                            "others_perc" = "KNON \n 2%",
#                            "r_perc" = "Others \n 9%"))

mean_data <- Simon_mass %>%
  mutate(Event_F = PMF_BA_full$color) %>%
  summarise(across(c(ii_perc, om_perc, ec_perc, gm_perc, ss_perc, others_perc, r_perc),\(x) mean(x, na.rm = TRUE))) %>%
  pivot_longer(cols = c(ii_perc, om_perc, ec_perc, gm_perc, ss_perc, others_perc, r_perc),
               names_to = "category", values_to = "value") %>%
  mutate(category = recode(category,
                           "ii_perc" = "II \n (10%)",
                           "om_perc" = "OM (57%)",
                           "ec_perc" = "EC (6%)",
                           "gm_perc" = "GM \n (14%)",
                           "ss_perc" = "SS (3%)",
                           "others_perc" = "KNON \n (2%)",
                           "r_perc" = "Others(9%)"))
# Create average pie chart 
mean_data <-mean_data %>% 
  arrange(desc(category)) %>%
  mutate(percentage = value ,
         ypos = cumsum(value) - 0.5* value)

# Crear gráfico de barras de PM2.5 ####
bar_plot <- ggplot(PMF_BA_full, aes(x = date2, y = `PM2,5`, fill = color)) +
  geom_bar(stat = "identity", width = 2.5) +
  scale_fill_manual(values = c("BB samples" = "darkred", "non-BB samples" = "darkblue")) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
  labs(x = "", y = "PM2.5 (µg/m³)", fill = "") +
  theme_minimal() +theme(legend.position = "top", 
                         plot.title = element_text(size = 20, face = "bold"),  # Tamaño del título
                         text = element_text(size = 24))

bar_plot



average_pie <- ggplot(data= mean_data, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  geom_label_repel(data = mean_data,
                   aes(label = category, x = .9, y=ypos), size = 8, nudge_x =.9,
                   box.padding = 0.1, point.padding = 0.5,
                   segment.color = 'grey50', show.legend = TRUE,
                   direction = "y") +  
  labs(title = "Chemical Profile") +
  theme_void() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 20, face = "bold"),  # Tamaño del título
        text = element_text(size = 20)) # Hide legend
average_pie

# Combinar gráficos: barras a la izquierda, tortas a la derecha
combined_plot <- (bar_plot | (average_pie )) + 
  plot_layout(widths = c(3.4, 2)) 

# Mostrar gráfico combinado
print(combined_plot)
n=4
png(filename="../Figures/Graphicalabstractreview2.png", res=300, height=n*5,width = n*13, units = "cm" )
print(combined_plot)
dev.off()

# polar maps####
library(openairmaps)
library(leaflet)
customIcon <- makeIcon(
  iconUrl = "/home/usuario/mdiaz/Documents/paper_facu/Figures/central.png",  
  iconWidth = 120,  
  iconHeight = 120
)

customIconships <- makeIcon(
  iconUrl = "/home/usuario/mdiaz/Documents/paper_facu/Figures/ship-flat-boat-by-Vexels.svg",  
  iconWidth = 120,  
  iconHeight = 120
)

customIcontree <- makeIcon(
  iconUrl = "/home/usuario/mdiaz/Documents/paper_facu/Figures/tree.png",  
  iconWidth = 120,  
  iconHeight = 120
)


customIconcars <- makeIcon(
  iconUrl = "/home/usuario/mdiaz/Documents/paper_facu/Figures/car-fleet-12792.png",  
  iconWidth = 120,  
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
minlat <- -60
maxlat <- 20
minlon <- -100
maxlon <- 0 
graphlimits_cwt  <-c(0, 60)
lonlatinc <-.5
# Generar el mapa de trayectorias ucustomIconcars# Generar el mapa de trayectorias usando trajMap
trajMap(selectByDate(subset(traj500, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat), 
                     start = "2019-08-26", end = "2019-08-30"),
        origin = TRUE,  
        grid.col = "transparent", 
        map.cols = "transparent",
        projection = "stereographic", 
        orientation = c(0, -65, 0), 
        parameters = NULL,  static = TRUE,
        static.nrow = TRUE,  control.position = "topright",
        control = NULL)

traj_plot<-trajPlot(
  selectByDate(
    subset(traj500, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat), 
    start = "2019-08-26", 
    end = "2019-08-30"
  ),
  plot.type="l",
  origin = FALSE,
  projection = "stereographic",
  orientation = c(0, -65, 0),  # Centrado en latitud -65
  parameters = NULL,
  
)
# Extraer el mapa base desde `traj_plot`
base_map <- traj_plot$map

# Agregar los datos de incendios como puntos al mapa
final_plot <- base_map +
  geom_point(data = fires_filtered, aes(x = longitude, y = latitude, color = datetime),
             size = 2, alpha = 0.7) +
  labs(color = "Fire Time") +  # Etiqueta para la leyenda de incendios
  theme_minimal()

###################33
library(openair)
library(ggplot2)

# Filtrado de incendios ya realizado
fires_filtered <- fires %>%
  filter(datetime >= start_filter & datetime <= end_filter)

# Configurar límites del mapa
minlon <- -80
maxlon <- -50
minlat <- -40
maxlat <- -20

# Crear el mapa de trayectorias usando trajLevel NO FUNCIONA
traj_data <- selectByDate(
  subset(traj500, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat), 
  start = "2019-08-26", 
  end = "2019-08-30"
)

traj_plot <- trajLevel(
  traj_data,
  lon = "lon",
  lat = "lat",
  pollutant = "height",
  grid.res = 1,
  map.cols = "lightgray",
  orientation = c(0, -65, 0),
  map.alpha = 0.7
)

# Convertir a ggplot para superposición
traj_ggplot <- traj_plot$plot

# Agregar los puntos de incendios al gráfico
final_plot <- traj_ggplot +
  geom_point(data = fires_filtered, aes(x = longitude, y = latitude, color = datetime),
             size = 2, alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", name = "Fire Time") +
  theme_minimal()

# Mostrar el gráfico
print(final_plot)

#########


# Mostrar el gráfico
print(final_plot)
library(sf)
library(ggplot2)

# Crear un objeto sf para las trayectorias
traj_sf <- st_as_sf(traj500, coords = c("lon", "lat"), crs = 4326) 

# Filtrar las trayectorias por latitud, longitud y rango de fechas
traj_filtered <- traj_sf %>%
  filter(
    lon >= minlon & lon <= maxlon,
    lat >= minlat & lat <= maxlat,
    date >= as.Date("2019-08-26") & date <= as.Date("2019-08-30")
  )

# Crear el mapa
ggplot() +
  geom_sf(data = traj_filtered, aes(color = factor(id)), size = 0.5) + # Color por ID de trayectoria
  scale_color_manual(values = c("darkred", "darkblue")) + # Colores personalizados
  coord_sf(xlim = c(minlon, maxlon), ylim = c(minlat, maxlat), expand = FALSE) + # Limites del mapa
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude", color = "Trajectory ID") +
  theme(
    panel.grid.major = element_line(color = "gray80"),  # Líneas de lat/lon
    panel.border = element_rect(color = "black", fill = NA), # Bordes del mapa
    legend.position = "top"
  )

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
  static = TRUE,
  static.nrow = NULL,
  progress = TRUE
)

# Dock sud -34.648961551867856, -58.34236555984067
# Central costanera -34.64608663021544, -58.344426375249924
# Central Puerto -34.57504109311285, -58.380487846322495

# #GA ####
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# library(ggrepel)
# 
# # Renaming categories in the DataFrame
# mean_data <- Simon_mass %>%
#   mutate(Event_F = PMF_BA_full$color) %>%
#   group_by(Event_F) %>%
#   summarise(across(c(inorganic_ions, organic_mass, elemental_C, geological_minerals, salt, others), mean, na.rm = TRUE)) %>%
#   pivot_longer(cols = c(inorganic_ions, organic_mass, elemental_C, geological_minerals, salt, others),
#                names_to = "category", values_to = "value") %>%
#   mutate(category = recode(category,
#                            "inorganic_ions" = "II",
#                            "organic_mass" = "OM",
#                            "elemental_C" = "EC",
#                            "geological_minerals" = "GM",
#                            "salt" = "SS",
#                            "others" = "others"))
# 
# # Define categories to be shown outside the pie chart
# 
# outside_labels <- c("EC","SS", "others")
# 
# 
# library(ggplot2)
# library(ggrepel)
# library(tidyverse)
# library(patchwork)
# # Crear gráfico de barras de PM2.5
# PMF_BA_full$color <- ifelse(PMF_BA_full$Event_F %in% c("SI", "SF","SO"), "BB samples", "non-BB samples")
# 
# PMF_BA_full$date=factor(PMF_BA_full$date)
# bar_plot <- ggplot(PMF_BA_full, aes(x = date, y = `PM2,5`, fill = color)) +
#   geom_bar(stat = "identity", width = 0.95) +
#   scale_fill_manual(values = c("BB samples" = "red", "non-BB samples" = "blue")) +
#   scale_x_discrete(breaks = levels(PMF_BA_full$date)[seq(1, length(levels(PMF_BA_full$date)), by = 4)]) +
#   labs(x = "Date", y = "PM2.5 (µg/m³)", fill = "Samples") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# 
# # # Create pie chart for events
# # event_data <-mean_data %>%
# #   filter(Event_F == "Event") %>%
# #   arrange(desc(category)) %>%
# #   mutate(percentage = value / sum(value) * 100,
# #          ypos = cumsum(value) - 0.5 * value) 
# # event_pie <- ggplot(event_data, aes(x = "", y = value, fill = category)) +
# #   geom_bar(stat = "identity", width = 1) +
# #   coord_polar("y") +
# #   geom_label_repel(data = event_data,
# #                    aes(label = category, x = 1.1), size = 4, nudge_x = 0.5,
# #                    box.padding = 0.5, point.padding = 0.5,
# #                    segment.color = 'grey50', show.legend = FALSE,
# #                    direction = "y") +  # Use direction "y" to ensure labels move along the y-axis
# #   labs(title = "BB samples") +
# #   theme_void() +
# #   theme(legend.position = "none") # Hide legend
# # event_pie
# library(ggplot2)
# library(dplyr)
# library(ggrepel)
# 
# # Preparar los datos para eventos
# event_data <- mean_data %>%
#   filter(Event_F == "Event") %>%
#   arrange(desc(category)) %>%
#   mutate(percentage = value / sum(value) * 100,
#          ypos = cumsum(value) - 0.5 * value)  # Calcular posiciones para las etiquetas
# 
# # Crear gráfico de torta con geom_label_repel()
# event_pie <- ggplot(event_data, aes(x = "", y = value, fill = category)) +
#   geom_bar(stat = "identity", width = 1) +
#   coord_polar("y") +
#   
#   # Colocar las etiquetas utilizando ypos y ajustes específicos para "others" y "SS"
#   geom_label_repel(aes(y = ypos, label = category), 
#                    size = 4, 
#                    nudge_x = ifelse(event_data$category %in% c("others", "SS"), 1, 0.5),  # Ajuste mayor para others y SS
#                    nudge_y = ifelse(event_data$category %in% c("others", "SS"), 0.2, 0),  # Ajustar verticalmente si es necesario
#                    box.padding = 0.5, 
#                    point.padding = 0.5,
#                    segment.color = 'grey50', 
#                    show.legend = FALSE,
#                    direction = "y") +  # Coloca las etiquetas a lo largo del eje Y
#   
#   labs(title = "BB samples") +
#   theme_void() +
#   theme(legend.position = "none")  # Ocultar leyenda
# 
# event_pie
# 
# 
# 
# 
# # Crear gráfico de pastel para eventos sin etiquetas
# no_event_pie <- ggplot(mean_data %>% filter(Event_F == "No Event"), aes(x = "", y = value, fill = category)) +
#   geom_bar(stat = "identity", width = 1) +
#   coord_polar("y") +
#   # Etiquetas dentro del gráfico para categorías que caben
#   geom_text(data = filter(mean_data %>% filter(Event_F == "No Event")),
#             aes(label = category), 
#             position = position_stack(vjust = 0.5), size = 4, 
#             hjust = 0.5) +
#   labs(title = "No Event") +
#   theme_void() +
#   theme(legend.position = "none") # Ocultar leyenda
# 
# # Combine plots: bars on the left, pies on the right
# combined_plot <- (bar_plot | (event_pie / no_event_pie))
# 
# # Show combined plot
# print(combined_plot)
# 
# 
# # Combinar gráficos: barras a la izquierda, tortas a la derecha
# combined_plot <- (bar_plot / (event_pie | no_event_pie))
# 
# # Mostrar gráfico combinado
# print(combined_plot)