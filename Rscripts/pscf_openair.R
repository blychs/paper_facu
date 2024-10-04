# Seteo de librerias ####
library(openair)
library(ggplot2)
library(lubridate)
library(readxl)
library(openair)
library(ggplot2)
library(lubridate)
library(readxl)
setwd("~/mdiaz//Documents/paper_facu/Rscripts")
source('opentrajutils.R')

# Seteo de parametros ####
# setwd("~/Documents/paper_facu/Rscripts/")
pathgraphs="Figures"
upper_windspeed=8.5
pscflims = c(0,1)
colorset="plasma"
pathtraj="../tdumps"
pattern ="tdump*"
bbpdatafile="../data_every_hour_obsv5.csv"

# Load datasets ####
bbpdata <- read.csv(bbpdatafile)
bbpdata$date <- as.POSIXct(bbpdata$date, tz='UTC')
bbpdata <- bbpdata[,-c(8:9,30:42)]
bbpdata$Sb_ng=bbpdata$Sb*1000
bbpdata$As_ng=bbpdata$As*1000
bbpdata$OC_EC = bbpdata$`C.Orgánico`/bbpdata$`C.Elemental`
bbpdata$KNON = (bbpdata$K - 0.6*bbpdata$Fe)
traj500<-read.tdumps(pathtraj = pathtraj, pattern=pattern)
traj500=traj500[traj500$receptor==3,]
trajconchem <- add.chem(bbpdata,traj500)
trajconchem = dplyr::rename(trajconchem, EC=C.Elemental)
trajconchem = dplyr::rename(trajconchem, OC=C.Orgánico)
trajconchem = dplyr::rename(trajconchem, TC=C.Total)
trajconchem = dplyr::rename(trajconchem, Na=Na.sol)

# define PCSF domain (Long range) #### 
minlat <- -60
maxlat <- 0
minlon <- -90
maxlon <- -25 
graphlimits_cwt  <-c(0, 60)
lonlatinc <- 0.5

# define PCSF domain (Long range) #### 
minlat <- -60
maxlat <- 20
minlon <- -90
maxlon <- 0 
graphlimits_cwt  <-c(0, 60)
lonlatinc <-.5


#72.8 DE TRAJ ES EL 75 DE LA MUESTRA <- chequear esto

# PSCF loop ####

keys = names(trajconchem)[c(14:42)]
print(keys)

for (title in keys) {
  print(title)
  png(paste0("../",pathgraphs,"/PSCF500m_75p_",title,".png"), width = 590 * 3, height =592* 3, res = 300)
  pscf<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                  pollutant = title , statistic = "pscf", limits = pscflims, percentile = 75,
                  projection = "stereographic",   orientation=c(0,-65,0), parameters = NULL,
                  cols = "heat",  border = NA, 
                  grid.cols = "grey20", auto.text =FALSE,  key.header=title,
                  origin = TRUE,
                  lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
  dev.off()
  }

title="EC"
pscfEC<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                pollutant = title , statistic = "pscf", limits = pscflims, percentile = 75,
                projection = "stereographic",   orientation=c(0,-65,0), parameters = NULL,
                cols = "heat",  border = NA, 
                grid.cols = "grey20", auto.text =FALSE, key.header =title, 
                origin = TRUE,
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
title="Fe"
pscfFe<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                  pollutant = title , statistic = "pscf", limits = pscflims, percentile = 75,
                  projection = "stereographic",  orientation=c(0,-65,0), parameters = NULL,
                  cols = "heat",  border = NA, 
                  grid.cols = "grey20", auto.text =FALSE, key.header =title, 
                  origin = TRUE,
                  lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
title="Na"
pscfNa<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                  pollutant = title , statistic = "pscf", limits = pscflims, percentile = 75,
                  projection = "stereographic",   orientation=c(0,-65,0), parameters = NULL,
                  cols = "heat",  border = NA, 
                  grid.cols = "grey20", auto.text =FALSE, key.header =title, 
                  origin = TRUE,
                  lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)

png(paste0("../",pathgraphs,"/pscfconjunto_",colorset,".png"), width = 717 * 5, height = 339* 5, res = 300)
print(pscfEC,split=c(1, 1, 3, 1))
print(pscfFe,split=c(2, 1, 3, 1), newpage=FALSE)
print(pscfNa,split=c(3, 1, 3, 1), newpage=FALSE)
dev.off()
dev.off()

library(patchwork)
selected_keys <- c("OC", "Fe", "Na") 
plots <- list()

for (i in seq_along(selected_keys)) {
  title <- selected_keys[i]
  
  # Aquí generas el gráfico de PSCF para cada key
  pscf <- trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                    pollutant = title, statistic = "pscf", limits = pscflims, percentile = 75,
                    projection = "stereographic", orientation = c(0, -65, 0), parameters = NULL,
                    cols = "heat", border = NA, 
                    grid.cols = "grey20", auto.text = FALSE, key.header = "", 
                    origin = TRUE,
                    lon.inc = lonlatinc, lat.inc = lonlatinc, min.bin = 2)
  
  # Convertir el gráfico a un objeto ggplot si es necesario
  # Aquí deberías ajustar según cómo se genere el gráfico (dependiendo de tu implementación)
  
  # Agregar anotaciones
  pscf <- pscf +
    annotate("text", x = Inf, y = Inf, label = letters[i], vjust = 1.5, hjust = 1.5, size = 5, fontface = "bold") +
    annotate("text", x = Inf, y = -Inf, label = title, vjust = -0.5, hjust = 1.5, size = 5, fontface = "italic")
  
  # Guardar el gráfico en la lista
  plots[[i]] <- pscf
}

# Combinar los gráficos en un solo panel usando patchwork
combined_plot <- wrap_plots(plots, ncol = 3)

# Guardar el gráfico combinado
png(paste0("../", pathgraphs, "/combined_pscf_plots.png"), width = 590 * 3, height = 592 * 3, res = 300)
print(combined_plot)
dev.off()

# PSCF EC ####
png(paste0(pathgraphs,"/PSCF500m_75p_EC.png"), width = 590 * 3, height =592* 3, res = 300)
pscf<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                pollutant = "EC", statistic = "pscf", 
                percentile = 75,
                projection = "stereographic",   orientation=c(0,-65,0), 
                parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
dev.off()

# PSCF EC Mel ####
png(paste0(pathgraphs,"/PSCF500m_75p_ECMel.png"), width = 590 * 3, height =592* 3, res = 300)
pscf<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                pollutant = "EC", statistic = "pscf", 
                percentile = 70,
                projection = "stereographic",   orientation=c(0,-65,0), limits=c(0.1,0.7),
                parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
dev.off()
# PSCF OC ####
png(paste0(pathgraphs,"/PSCF500m_75p_OC.png"), width = 590 * 3, height =592* 3, res = 300)
pscf<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                pollutant = "OC", statistic = "pscf", 
                percentile = 75,
                projection = "stereographic",   orientation=c(0,-65,0), 
                parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
dev.off()

# PSCF K ####
png(paste0(pathgraphs,"/PSCF500m_75p_K.png"), width = 590 * 3, height =592* 3, res = 300)
pscf<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                pollutant = "K", statistic = "pscf", limits = pscflims,
                percentile = 75,
                projection = "stereographic",   orientation=c(0,-65,0), 
                parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
dev.off()

# PSCF KNON ####
png(paste0(pathgraphs,"/PSCF500m_75p_KNON.png"), width = 590 * 3, height =592* 3, res = 300)
pscf<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                pollutant = "KNON", statistic = "pscf", limits = pscflims,
                percentile = 75,
                projection = "stereographic",   orientation=c(0,-65,0), 
                parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
dev.off()

# PSCF Cl ####
png(paste0(pathgraphs,"/PSCF500m_75p_Cl.png"), width = 590 * 3, height =592* 3, res = 300)
pscf<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                pollutant = "Cl", statistic = "pscf", 
                percentile = 75,
                projection = "stereographic",   orientation=c(0,-65,0), 
                parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
dev.off()

# PSCF Na sol ####
png(paste0(pathgraphs,"/PSCF500m_75p_Nasol.png"), width = 590 * 3, height =592* 3, res = 300)
pscf<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                pollutant = "Na.sol", statistic = "pscf", 
                percentile = 75,
                projection = "stereographic",   orientation=c(0,-65,0), 
                parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "black", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "grey60", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
dev.off()

# PSCF Na sol ####
png(paste0(pathgraphs,"/PSCF500m_75p_Nasol.png"), width = 590 * 3, height =592* 3, res = 300)
pscf<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                pollutant = "Na.sol", statistic = "pscf", 
                percentile = 75,
                projection = "stereographic",   orientation=c(0,-65,0), parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
dev.off()

# PSCF Cu ####
png(paste0(pathgraphs,"/PSCF500m_75p_Cu.png"), width = 590 * 3, height =592* 3, res = 300)
pscf<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                pollutant = "Cu", statistic = "pscf", 
                percentile = 75,
                projection = "stereographic",   orientation=c(0,-65,0), 
                parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
dev.off()

# PSCF PM2.5 ####
png(paste0(pathgraphs,"/PSCF500m_75p_PM25.png"), width = 590 * 3, height =592* 3, res = 300)
pscf<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                pollutant = "PM2.5", statistic = "pscf", 
                percentile = 75,
                projection = "stereographic",   orientation=c(0,-65,0), 
                parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
dev.off()

# PSCF Ti ####
png(paste0(pathgraphs,"/PSCF500m_75p_Ti.png"), width = 590 * 3, height =592* 3, res = 300)
pscf<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                pollutant = "Ti", statistic = "pscf", 
                percentile = 75,
                projection = "stereographic",   orientation=c(0,-65,0), 
                parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
dev.off()

# PSCF Fe ####
png(paste0(pathgraphs,"/PSCF500m_75p_Fe.png"), width = 590 * 3, height =592* 3, res = 300)
pscf<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                pollutant = "Fe", statistic = "pscf", 
                percentile = 75,
                projection = "stereographic",   orientation=c(0,-65,0), 
                parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
dev.off()


png(paste0(pathgraphs,"/PSCF500m_PM25eventos_",title,".png"), width = 590 * 3, height =592* 3, res = 300)
pscf<-trajLevel(selectByDate(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat), 
                             start = "2019-08-01", end = "2019-08-31"),
                pollutant = title , statistic = "pscf", limits = pscflims, percentile = 60,
                projection = "stereographic",   orientation=c(0,-65,0), parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
dev.off()

pscf<-trajLevel(selectByDate(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat), 
                             start = "2019-08-01", end = "2019-08-31"),
                pollutant = title , statistic = "cwt",
                projection = "stereographic",   orientation=c(0,-65,0), parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
pscf<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                pollutant = title , statistic = "pscf", limits = pscflims, percentile = 70,
                projection = "stereographic",   orientation=c(0,-65,0), parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)

pscf<-trajLevel(selectByDate(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat), 
                             start = "2019-07-22", end = "2019-09-29"),
                pollutant = title , statistic = "cwt",  
                projection = "stereographic",   orientation=c(0,-65,0), parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
pscf<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                pollutant = title , statistic = "cwt",  
                projection = "stereographic",   orientation=c(0,-65,0), parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
# 
pscf<-trajLevel(selectByDate(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat), 
                            start = "2019-08-07", end = "2019-08-28"),
               pollutant = "OC" , statistic = "pscf", limits = pscflims, percentile = 70,
               projection = "stereographic",   orientation=c(0,-65,0), parameters = NULL,
               cols = "heat", smooth = F,  border = NA, 
               grid.cols = "grey40", auto.text =FALSE, key.header ="", 
               origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
               lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)

pscf<-trajLevel(selectByDate(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat), 
                             start = "2019-09-06", end = "2019-09-27"),
                pollutant = "OC" , statistic = "pscf", limits = pscflims, percentile = 70,
                projection = "stereographic",   orientation=c(0,-65,0), parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)

pscf<-trajPlot(selectByDate(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat), 
                             start = "2019-04-02", end = "2019-04-10"),
               origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
               projection = "stereographic",   orientation=c(0,-65,0), parameters = NULL)


pscf<-trajPlot(selectByDate(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat), 
                            start = "2019-08-12", end = "2019-08-30"),
               origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
               projection = "stereographic",   orientation=c(0,-65,0), parameters = NULL)

trajPlot(selectByDate(subset(traj500, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat), 
                      start = "2019-08-25", end = "2019-09-01"),
         origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
         projection = "mercator",   orientation=c(0,-65,0), parameters = NULL)
pscf<-trajLevel(selectByDate(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat), 
                             start = "2020-01-01", end = "2020-01-23"),
                pollutant = "OC" , statistic = "pscf", limits = pscflims, percentile = 60,
                projection = "stereographic",   orientation=c(0,-65,0), parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
# cluster ####
clust <- trajCluster(traj500, method = "Euclid", n.cluster= 6, col = "Set2", origin = FALSE,
                     projection = "stereographic",   orientation=c(0,-65,0), parameters = NULL,
                     map.cols = openColours("Paired", 10))

clust <- trajCluster(traj500, method = "Angle", n.cluster= 6, col = "Set2", origin = FALSE,
                     projection = "stereographic",   orientation=c(0,-65,0), parameters = NULL,
                     map.cols = openColours("Paired", 10), 
                     type = "season",hemisphere="southern")
clustdata=subset(clust$data, hour.inc == 0)

datamergeclust=merge(PMF_BA_full,, by = "date")
trendLevel(datamergeclust, pollutant = "v2.5", type = "cluster", layout = c(6, 1),
           cols = "increment")
