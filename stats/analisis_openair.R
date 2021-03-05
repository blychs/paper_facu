library(readxl)
library(ggplot2)
library(openair)
library(tidyverse)

datos <- read_excel('Documents/openair_var.xlsx', sheet = 'calculos intermedios')
datos_meteo <- read.csv('doctorado/datos_smn/meteo_observatorio_2019-2020_calmacomoNA_paraR.csv')
View(datos_meteo)
View(datos)

day_array <- c()
for (i in 1:120){
  day_array <- c(day_array, (i*24 + 1)) 
}

print(day_array)
for (i in day_array){
  windRose(datos_meteo[i: (i+23),], key.header = (datos_meteo['Fecha'][i,]))
  }
