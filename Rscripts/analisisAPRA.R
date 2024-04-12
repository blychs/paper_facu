## Datos de las estaciones de monitoreo del Gobierno de la Ciudad

### Unidades

* NO, NO<sub>2</sub>, NO<sub>x</sub> y SO<sub>2</sub> en ppb \
* CO en ppm \
* PM<sub>10</sub> en μg/m^3^ \

### Formato de fecha 

YYYY-MM-DD HH:mm:ss

### Datos faltantes

* NA para s/d \

Los valores menores al limite de deteccion (LD) fueron reemplazados por un numero random menor a LD \

* LD_CO=0.05 \
* LD_NO<sub>x</sub>=4 \

## Observaciones
* Parque Centenario no mide SO<sub>2</sub> asi que fueron reemplazados por NA

### Leer y Filtrar datos de una única estación


library(readr)
library(lubridate)
workingdirectory = "C:/Users/User/Desktop/datoscalidadcnea"
setwd(workingdirectory)
datadirectory ="DatosAPRA"    
source("codes/module_graficos.R")
pathgraphs="C:/Users/User/Desktop/datoscalidadcnea/graficos/APRA/Boxplot"

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

png(paste0(pathgraphs,"CO_Boxplot_Centenario.png"), width = 8 * 300, height = 5* 300, res = 300)
boxplotbyMonth(dataestacion, pollutant = "CO")
dev.off()

png(paste0(pathgraphs,"NOX_Boxplot_Centenario.png"), width = 8 * 300, height = 5* 300, res = 300)
boxplotbyMonth(dataestacion, pollutant = "NOX")
dev.off()

png(paste0(pathgraphs,"NO_Boxplot_Centenario.png"), width = 8 * 300, height = 5* 300, res = 300)
boxplotbyMonth(dataestacion, pollutant = "NO")
dev.off()

png(paste0(pathgraphs,"NO2_Boxplot_Centenario.png"), width = 8 * 300, height = 5* 300, res = 300)
boxplotbyMonth(dataestacion, pollutant = "NO2")
dev.off()

png(paste0(pathgraphs,"PM10_Boxplot_Centenario.png"), width = 8 * 300, height = 5* 300, res = 300)
boxplotbyMonth(dataestacion, pollutant = "PM10")
dev.off()

#Estacion La Boca
dataestacion = datosAPRA[datosAPRA$estacion=="La Boca",]

png(paste0(pathgraphs,"CO_Boxplot_La Boca.png"), width = 8 * 300, height = 5* 300, res = 300)
boxplotbyMonth(dataestacion, pollutant = "CO")
dev.off()

png(paste0(pathgraphs,"NOX_Boxplot_La Boca.png"), width = 8 * 300, height = 5* 300, res = 300)
boxplotbyMonth(dataestacion, pollutant = "NOX")
dev.off()

png(paste0(pathgraphs,"NO_Boxplot_La Boca.png"), width = 8 * 300, height = 5* 300, res = 300)
boxplotbyMonth(dataestacion, pollutant = "NO")
dev.off()

png(paste0(pathgraphs,"NO2_Boxplot_La Boca.png"), width = 8 * 300, height = 5* 300, res = 300)
boxplotbyMonth(dataestacion, pollutant = "NO2")
dev.off()

png(paste0(pathgraphs,"PM10_Boxplot_La Boca.png"), width = 8 * 300, height = 5* 300, res = 300)
boxplotbyMonth(dataestacion, pollutant = "PM10")
dev.off()

png(paste0(pathgraphs,"So2_Boxplot_La Boca.png"), width = 8 * 300, height = 5* 300, res = 300)
boxplotbyMonth(dataestacion, pollutant = "SO2")
dev.off()

#Estacion Cordoba
dataestacion = datosAPRA[datosAPRA$estacion=="Cordoba",]

png(paste0(pathgraphs,"CO_Boxplot_Cordoba.png"), width = 8 * 300, height = 5* 300, res = 300)
boxplotbyMonth(dataestacion, pollutant = "CO")
dev.off()

png(paste0(pathgraphs,"NOX_Boxplot_Cordoba.png"), width = 8 * 300, height = 5* 300, res = 300)
boxplotbyMonth(dataestacion, pollutant = "NOX")
dev.off()

png(paste0(pathgraphs,"NO_Boxplot_Cordoba.png"), width = 8 * 300, height = 5* 300, res = 300)
boxplotbyMonth(dataestacion, pollutant = "NO")
dev.off()

png(paste0(pathgraphs,"NO2_Boxplot_Cordoba.png"), width = 8 * 300, height = 5* 300, res = 300)
boxplotbyMonth(dataestacion, pollutant = "NO2")
dev.off()

png(paste0(pathgraphs,"PM10_Boxplot_Cordoba.png"), width = 8 * 300, height = 5* 300, res = 300)
boxplotbyMonth(dataestacion, pollutant = "PM10")
dev.off()

png(paste0(pathgraphs,"So2_Boxplot_Cordoba.png"), width = 8 * 300, height = 5* 300, res = 300)
boxplotbyMonth(dataestacion, pollutant = "SO2")
dev.off()

### Para todos los calculos de medias y medianas!!#

#cargar datos de gases (23-02 al 25-05)

# usamos datosAPRA
#Estacion Cordoba
d_APRA = datosAPRA[datosAPRA$estacion=="Cordoba",]

# Historico feb 2019 a mayo 20
# Con la funcion selectbyDate fijas los periodos (Dejo uno hecho como ejemplo)
dathist=LD
LD=selectByDate(d_APRA, start = "2019-02-23", end = "2020-03-15")
# Con summary te tira todos los valores medios y medianas, copiar cada valor en el excel online
summary(LD)
#Con timeAverage calculas los valores medios diarios (Dejo un periodo de ejemplo, hacer el resto)
LDdaily=timeAverage(LD,avg.time = "day", statistic = "mean")
# Idem con hourly pero con daily: Con summary te tira todos los valores medios y medianas
summary(LDdaily)


## MAM (2019) 20-03 al 25-05
# Con la funcion selectbyDate fijas los periodos (Dejo uno hecho como ejemplo)
MAM=LD
LD=selectByDate(d_APRA, start = "2019-03-20", end = "2019-05-25")
# Con summary te tira todos los valores medios y medianas, copiar cada valor en el excel online
summary(LD)
#Con timeAverage calculas los valores medios diarios (Dejo un periodo de ejemplo, hacer el resto)
LDdaily=timeAverage(LD,avg.time = "day", statistic = "mean")
# Idem con hourly pero con daily: Con summary te tira todos los valores medios y medianas
summary(LDdaily)


## BLD  (2020) 16-02 al 16-03
# Con la funcion selectbyDate fijas los periodos (Dejo uno hecho como ejemplo)
BLD=LD
LD=selectByDate(d_APRA, start = "2020-02-16", end = "2020-03-16")
# Con summary te tira todos los valores medios y medianas, copiar cada valor en el excel online
summary(LD)
#Con timeAverage calculas los valores medios diarios (Dejo un periodo de ejemplo, hacer el resto)
LDdaily=timeAverage(LD,avg.time = "day", statistic = "mean")
summary(LDdaily)



## LD  (2020) 20-03 al 13-04
# Con la funcion selectbyDate fijas los periodos (Dejo uno hecho como ejemplo)
LD=LD
LD=selectByDate(d_APRA, start = "2020-03-20", end = "2020-04-13")
# Con summary te tira todos los valores medios y medianas, copiar cada valor en el excel online
summary(LD)
#Con timeAverage calculas los valores medios diarios (Dejo un periodo de ejemplo, hacer el resto)
LDdaily=timeAverage(LD,avg.time = "day", statistic = "mean")
# Idem con hourly pero con daily: Con summary te tira todos los valores medios y medianas
summary(LDdaily)


## PLD  (2020) 14-04 al 25-05
# Con la funcion selectbyDate fijas los periodos (Dejo uno hecho como ejemplo)
PLD=LD
LD=selectByDate(d_APRA, start = "2020-04-14", end = "2020-05-25")
# Con summary te tira todos los valores medios y medianas, copiar cada valor en el excel online
summary(LD)
#Con timeAverage calculas los valores medios diarios (Dejo un periodo de ejemplo, hacer el resto)
LDdaily=timeAverage(LD,avg.time = "day", statistic = "mean")
# Idem con hourly pero con daily: Con summary te tira todos los valores medios y medianas
summary(LDdaily)


#Estacion La Boca
d_APRA = datosAPRA[datosAPRA$estacion=="La Boca",]

# Historico feb 2019 a mayo 20
# Con la funcion selectbyDate fijas los periodos (Dejo uno hecho como ejemplo)
dathist=LD
LD=selectByDate(d_APRA, start = "2019-02-23", end = "2020-03-15")
# Con summary te tira todos los valores medios y medianas, copiar cada valor en el excel online
summary(LD)
#Con timeAverage calculas los valores medios diarios (Dejo un periodo de ejemplo, hacer el resto)
LDdaily=timeAverage(LD,avg.time = "day", statistic = "mean")
# Idem con hourly pero con daily: Con summary te tira todos los valores medios y medianas
summary(LDdaily)

## MAM (2019) 20-03 al 25-05
# Con la funcion selectbyDate fijas los periodos (Dejo uno hecho como ejemplo)
MAM=LD
LD=selectByDate(d_APRA, start = "2019-03-20", end = "2019-05-25")
# Con summary te tira todos los valores medios y medianas, copiar cada valor en el excel online
summary(LD)
#Con timeAverage calculas los valores medios diarios (Dejo un periodo de ejemplo, hacer el resto)
LDdaily=timeAverage(LD,avg.time = "day", statistic = "mean")
# Idem con hourly pero con daily: Con summary te tira todos los valores medios y medianas
summary(LDdaily)


## BLD  (2020) 16-02 al 16-03
# Con la funcion selectbyDate fijas los periodos (Dejo uno hecho como ejemplo)
BLD=LD
LD=selectByDate(d_APRA, start = "2020-02-16", end = "2020-03-16")
# Con summary te tira todos los valores medios y medianas, copiar cada valor en el excel online
summary(LD)
#Con timeAverage calculas los valores medios diarios (Dejo un periodo de ejemplo, hacer el resto)
LDdaily=timeAverage(LD,avg.time = "day", statistic = "mean")
# Idem con hourly pero con daily: Con summary te tira todos los valores medios y medianas
summary(LDdaily)
#LDhourly=timeAverage(LD,avg.time = "hour", statistic = "mean")
#summary(LDhourly)


## LD  (2020) 20-03 al 13-04
# Con la funcion selectbyDate fijas los periodos (Dejo uno hecho como ejemplo)
LD=LD
LD=selectByDate(d_APRA, start = "2020-03-20", end = "2020-04-13")
# Con summary te tira todos los valores medios y medianas, copiar cada valor en el excel online
summary(LD)
#Con timeAverage calculas los valores medios diarios (Dejo un periodo de ejemplo, hacer el resto)
LDdaily=timeAverage(LD,avg.time = "day", statistic = "mean")
# Idem con hourly pero con daily: Con summary te tira todos los valores medios y medianas
summary(LDdaily)


## PLD  (2020) 14-04 al 25-05
# Con la funcion selectbyDate fijas los periodos (Dejo uno hecho como ejemplo)
PLD=LD
LD=selectByDate(d_APRA, start = "2020-04-14", end = "2020-05-25")
# Con summary te tira todos los valores medios y medianas, copiar cada valor en el excel online
summary(LD)
#Con timeAverage calculas los valores medios diarios (Dejo un periodo de ejemplo, hacer el resto)
LDdaily=timeAverage(LD,avg.time = "day", statistic = "mean")
# Idem con hourly pero con daily: Con summary te tira todos los valores medios y medianas
summary(LDdaily)

#Estacion Centenario
d_APRA = datosAPRA[datosAPRA$estacion=="Centenario",]
#View(d_APRA)

# Historico feb 2019 a mayo 20
# Con la funcion selectbyDate fijas los periodos (Dejo uno hecho como ejemplo)
dathist=LD
LD=selectByDate(d_APRA, start = "2019-02-23", end = "2020-03-15")
# Con summary te tira todos los valores medios y medianas, copiar cada valor en el excel online
summary(LD)
#Con timeAverage calculas los valores medios diarios (Dejo un periodo de ejemplo, hacer el resto)
LDdaily=timeAverage(LD,avg.time = "day", statistic = "mean")
# Idem con hourly pero con daily: Con summary te tira todos los valores medios y medianas
summary(LDdaily)
#me falta agregar la parte de pm2.5

## MAM (2019) 20-03 al 25-05
# Con la funcion selectbyDate fijas los periodos (Dejo uno hecho como ejemplo)
MAM=LD
LD=selectByDate(d_APRA, start = "2019-03-20", end = "2019-05-25")
# Con summary te tira todos los valores medios y medianas, copiar cada valor en el excel online
summary(LD)
#Con timeAverage calculas los valores medios diarios (Dejo un periodo de ejemplo, hacer el resto)
LDdaily=timeAverage(LD,avg.time = "day", statistic = "mean")
# Idem con hourly pero con daily: Con summary te tira todos los valores medios y medianas
summary(LDdaily)


## BLD  (2020) 16-02 al 16-03
# Con la funcion selectbyDate fijas los periodos (Dejo uno hecho como ejemplo)
BLD=LD
LD=selectByDate(d_APRA, start = "2020-02-16", end = "2020-03-16")
# Con summary te tira todos los valores medios y medianas, copiar cada valor en el excel online
summary(LD)
#Con timeAverage calculas los valores medios diarios (Dejo un periodo de ejemplo, hacer el resto)
LDdaily=timeAverage(LD,avg.time = "day", statistic = "mean")
# Idem con hourly pero con daily: Con summary te tira todos los valores medios y medianas
summary(LDdaily)
#LDhourly=timeAverage(LD,avg.time = "hour", statistic = "mean")
#summary(LDhourly)


## LD  (2020) 20-03 al 13-04
# Con la funcion selectbyDate fijas los periodos (Dejo uno hecho como ejemplo)
LD=LD
LD=selectByDate(d_APRA, start = "2020-03-20", end = "2020-04-13")
# Con summary te tira todos los valores medios y medianas, copiar cada valor en el excel online
summary(LD)
#Con timeAverage calculas los valores medios diarios (Dejo un periodo de ejemplo, hacer el resto)
LDdaily=timeAverage(LD,avg.time = "day", statistic = "mean")
# Idem con hourly pero con daily: Con summary te tira todos los valores medios y medianas
summary(LDdaily)


## PLD  (2020) 14-04 al 25-05
# Con la funcion selectbyDate fijas los periodos (Dejo uno hecho como ejemplo)
PLD=LD
LD=selectByDate(d_APRA, start = "2020-04-14", end = "2020-05-25")
# Con summary te tira todos los valores medios y medianas, copiar cada valor en el excel online
summary(LD)
#Con timeAverage calculas los valores medios diarios (Dejo un periodo de ejemplo, hacer el resto)
LDdaily=timeAverage(LD,avg.time = "day", statistic = "mean")
# Idem con hourly pero con daily: Con summary te tira todos los valores medios y medianas
summary(LDdaily)

