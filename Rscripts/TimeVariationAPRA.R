


library(readr)
library(lubridate)


workingdirectory = "C:/Users/User/Desktop/datoscalidadcnea"
setwd(workingdirectory)
datadirectory ="DatosAPRA"    
source("codes/module_graficos.R")
pathgraphs="C:/Users/User/Desktop/datoscalidadcnea/graficos/APRA/Timevariation/"

datosAPRA <- read_csv("DatosAPRA/datosAPRAlimpios.csv", 
                      col_types = cols(SO2 = col_double(), 
                                       date = col_datetime(format = "%Y-%m-%d %H:%M:%S")))
# normalizamos las unidades 
datosAPRA$NO=datosAPRA$NO/1000 # ppb a ppm
datosAPRA$NO2=datosAPRA$NO2/1000 # ppb a ppm
datosAPRA$NOX=datosAPRA$NOX/1000 # ppb a ppm
datosAPRA$SO2=datosAPRA$SO2/1000 # ppb a ppm

#View(datosAPRA)

#Estacion Centenario
dataestacion = datosAPRA[datosAPRA$estacion=="Centenario",]

png(paste0(pathgraphs,"CO_Timevariation_Centenario.png"), width = 8 * 300, height = 5* 300, res = 300)
timeVariation(dataestacion, pollutant = "CO", statistic = "median")
dev.off()

png(paste0(pathgraphs,"NOX_Timevariation_Centenario.png"), width = 8 * 300, height = 5* 300, res = 300)
timeVariation(dataestacion, pollutant = "NOX", statistic = "median")
dev.off()

png(paste0(pathgraphs,"NO_Timevariation_Centenario.png"), width = 8 * 300, height = 5* 300, res = 300)
timeVariation(dataestacion, pollutant = "NO", statistic = "median")
dev.off()

png(paste0(pathgraphs,"NO2_Timevariation_Centenario.png"), width = 8 * 300, height = 5* 300, res = 300)
timeVariation(dataestacion, pollutant = "NO2", statistic = "median")
dev.off()

png(paste0(pathgraphs,"PM10_Timevariation_Centenario.png"), width = 8 * 300, height = 5* 300, res = 300)
timeVariation(dataestacion, pollutant = "PM10", statistic = "median")
dev.off()



