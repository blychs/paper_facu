library(openair)
library(ggplot2)
library(lubridate)
library(readxl)
source('Rscripts/opentrajutils.R')

#setwd("/home/pablo/doctorado/paper_facu")
pathtraj="tdumps"

pathplots=pathgraphs

pattern ="tdump*"
traj500<-read.tdumps(pathtraj = pathtraj, pattern=pattern)
traj500=traj500[traj500$receptor==3,]

trajconchem <- add.chem(mergeddata,traj500) 

#Long range
minlat <- -60
maxlat <- 0
minlon <- -90
maxlon <- -25 
graphlimits_cwt  <-c(0, 60)

lonlatinc <- 0.5
#72.8 DE TRAJ ES EL 75 DE LA MUESTRA
png(paste0(pathgraphs,"/PSCF500m_75p_EC.png"), width = 590 * 3, height =592* 3, res = 300)
pscf<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                pollutant = "C.Elemental", statistic = "pscf", 
                percentile = 75,
                projection = "stereographic",   orientation=c(0,-65,0), 
                parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
dev.off()

png(paste0(pathgraphs,"/PSCF500m_75p_K.png"), width = 590 * 3, height =592* 3, res = 300)
pscf<-trajLevel(subset(trajconchem, lon >= minlon & lon <= maxlon & lat >= minlat & lat <= maxlat),
                pollutant = "K", statistic = "pscf", 
                percentile = 75,
                projection = "stereographic",   orientation=c(0,-65,0), 
                parameters = NULL,
                cols = "heat", smooth = F,  border = NA, 
                grid.cols = "grey40", auto.text =FALSE, key.header ="", 
                origin = TRUE,  grid.col = "transparent", map.cols = "transparent",
                lon.inc = lonlatinc , lat.inc = lonlatinc , min.bin = 2)
dev.off()
