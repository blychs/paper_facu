# Seteo de librerias ####
library(openair)
library(ggplot2)
library(lubridate)
library(readxl)
library(openair)
library(ggplot2)
library(lubridate)
library(readxl)
source('Rscripts/opentrajutils.R')

# Seteo de parametros ####
#setwd("/home/pablo/doctorado/paper_facu")
pathgraphs="Figures"
upper_windspeed=8.5
pscflims = c(0,1)
colorset="plasma"
pathtraj="tdumps"
pattern ="tdump*"
bbpdatafile="data/data_every_hour_obs_eventM.csv"

# Load datasets ####
bbpdata <- read.csv(bbpdatafile)
bbpdata$date <- as.POSIXct(bbpdata$date, tz='UTC')
bbpdata <- bbpdata[,-c(8,31:42)]
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

# define PCSF domain (Long range) #### 
minlat <- -60
maxlat <- 0
minlon <- -90
maxlon <- -25 
graphlimits_cwt  <-c(0, 60)
lonlatinc <- 0.5

#72.8 DE TRAJ ES EL 75 DE LA MUESTRA <- chequear esto

# PSCF loop ####
keys = names(trajconchem)[c(14:44, 46,52)]
print(keys)
for (title in keys) {
  print(title)
  png(paste0(pathgraphs,"/PSCF500m_75p_",title,".png"), width = 590 * 3, height =592* 3, res = 300)
  pscf<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                  pollutant = title , statistic = "pscf", limits = pscflims, percentile = 75,
                  projection = "stereographic",   orientation=c(0,-65,0), parameters = NULL,
                  cols = "heat", smooth = F,  border = NA, 
                  grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                  origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                  lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
  dev.off()
  }


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