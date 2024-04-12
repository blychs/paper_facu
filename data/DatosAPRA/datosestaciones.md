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

~~~
library(readr)
source("../codes/module_graficos.R")

datosAPRA <- read_csv("datosAPRAlimpios.csv", 
    col_types = cols(SO2 = col_double(), 
        date = col_datetime(format = "%Y-%m-%d %H:%M:%S")))
dataestacion = datosAPRA[datosAPRA$estacion=="Cordoba",]
boxplotbyMonth(dataestacion, pollutant = "CO")
dataestacion = datosAPRA[datosAPRA$estacion=="La Boca",]
boxplotbyMonth(dataestacion, pollutant = "CO")
~~~

