library(readxl)
library(openair)
library(readr) 
library(ranger)
library (openair)
library(plyr)
library(dplyr)
library(rmweather)
library(ggplot2)
library(Dict)
source("Rscripts/module_metales.R")
# library(doParallel)
# registerDoParallel(cores=6)
# extrafont::font_import()
extrafont::loadfonts()
# Cargar datos ####
BA_events_testM <- read_excel("BA_events_testM.xlsx",
                              col_types = c("date", "skip","numeric", "numeric", "numeric",
                                            "skip", "text","numeric"))
BA_events_testM$date <- as.POSIXct(BA_events_testM$date, tz='UTC')
BA_events_testM$Event_M <- as.factor(BA_events_testM$Event_M)
BA_events_testM$Lote <- as.factor(BA_events_testM$Lote)
BA_events_testM <- BA_events_testM[,c(-6)]#tiro el lote porque ya lo tengo 

volumen_tish <- read_excel("data/volumen_tish_corregido2_espero_corregido3.xlsx")
LD <- read_excel("data/Limites de detección.xlsx")

# Planilla_conjunta_Metales <- read_excel("data/Planilla conjunta Metales.xlsx", sheet = "Sheet1", na = c("n.d", "d"))
# Planilla_conjunta_Metales$Lote=as.factor(Planilla_conjunta_Metales$Lote)
# Planilla_conjunta_Metales$Blanco=as.factor(Planilla_conjunta_Metales$Blanco)

Planilla_conjunta_Metales <- read_excel("data/Planilla conjunta Metales.xlsx", sheet = "Sheet1")
# completando la serie temporal con los valores de limites de detección
# Planilla_conjunta_Metales <- Planilla_conjunta_Metales %>%  mutate_at(vars(5:28), ~ replace(., . =="n.d", 0))
nd=Planilla_conjunta_Metales$Be[1]
Planilla_conjunta_Metales <- Planilla_conjunta_Metales %>% mutate(across(.cols = 5:28, .fns = ~ {
  compuesto <- as.character(colnames(.))
  limites_compuesto <- LD$LD[match(compuesto, LD$Compuesto)]
  replace(., . == nd, limites_compuesto)
}))
Planilla_conjunta_Metales$Lote=as.factor(Planilla_conjunta_Metales$Lote)
Planilla_conjunta_Metales$Blanco=as.factor(Planilla_conjunta_Metales$Blanco)

# summary(Planilla_conjunta_Metales)

# Procesamiento de datos de ug/l a ug ####
ugtotales=merge(Planilla_conjunta_Metales, volumen_tish, by="ID",all.x = TRUE)
ugtotales[,c(5:28)]=ugtotales[,c(5:28)]*ugtotales$volumen_aforo/1000*ugtotales$area_material/ugtotales$area_analizada

# Separo cuales son filtros y cuales son Blancos porque los filtros tienen asociada una fecha y un volumen
Filtros=ugtotales[which(ugtotales$Blanco==0),]
Filtros=merge(Filtros,BA_events_testM,by="ID")

Blancos=ugtotales[which(ugtotales$Blanco==1),]


# Calculo de Blancos #### 
selecciondeblanco="promediosinminmax" # promedio, promedioporlote, promediosinmaxnimin
Blanco_para_restar = calculo_de_blanco(Blancos, selecciondeblanco)

# Resto el blanco y divido por el volumen
for (i in seq(dim(Filtros)[1])){ 
  Lote=as.numeric(as.character(Filtros$Lote[i]))
  Filtros[i,c(5:28)] = (Filtros[i,c(5:28)] - Blanco_para_restar[which(as.numeric(as.character(Blanco_para_restar$Lote))==Lote),-25])/Filtros$`Vol (m3)`[i]
}
# Acá está reemplazando por 0 y lo estamos reemplazando por un valor dentro del limite de detección en PMF_BA_full
Filtros <- Filtros %>%  mutate_at(vars(5:28), ~ replace(., . < 0, 0))

ggplot(Filtros)+geom_point(aes(x=date, y=Ca, color="Ca"))+geom_line(aes(x=date, y=Ca, color="Ca"))+
  geom_point(aes(x=date, y=Ca, color="Ca PMFfull"),data=PMF_BA_full)+geom_line(aes(x=date, y=Ca, color="Ca PMFfull"),
                                                                             data=PMF_BA_full)+
  geom_point(aes(x=date, y=K, color="Kblancop"))+geom_line(aes(x=date, y=K, color="Kblancop")) +
  geom_point(aes(x=date, y=K, color="K PMFfull"),data=PMF_BA_full)+geom_line(aes(x=date, y=K, color="K PMFfull"),
                                                                     data=PMF_BA_full)
ggplot(Filtros)+
  geom_point(aes(x=date, y=K, color="Kblancop"))+geom_line(aes(x=date, y=K, color="Kblancop")) +
  geom_point(aes(x=date, y=K, color="K PMFfull"),data=PMF_BA_full)+geom_line(aes(x=date, y=K, color="K PMFfull"),
                                                                             data=PMF_BA_full)

ggplot(Filtros)+geom_point(aes(x=date, y=Ca, color="Ca"))+geom_line(aes(x=date, y=Ca, color="Ca"))+
  geom_point(aes(x=date, y=Ca, color="Ca PMFfull"),data=PMF_BA_full)+geom_line(aes(x=date, y=Ca, color="Ca PMFfull"),
                                                                               data=PMF_BA_full)

ggplot(Filtros)+geom_point(aes(x=date, y=Ca, color="Ca"))+geom_line(aes(x=date, y=Ca, color="Ca"))+
  geom_point(aes(x=date, y=K, color="Kblancop"))+geom_line(aes(x=date, y=K, color="Kblancop"))

ggplot(Filtros)+geom_point(aes(x=date, y=Ca, color="Ca"))+geom_line(aes(x=date, y=Ca, color="Ca"))+
  geom_point(aes(x=date, y=K, color="K"),data=PMF_BA_full)+geom_line(aes(x=date, y=K, color="K"), 
                                                                     data=PMF_BA_full)

ggplot(Filtros)+geom_point(aes(x=date, y=Ca, color="Ca"))+geom_line(aes(x=date, y=Ca, color="Ca"))+
  geom_point(aes(x=date, y=K, color="Kblancop"))+geom_line(data=PMF_BA_full,aes(x=date, y=`C Orgánico`/10, color="OC"))


ggplot(Filtros)+geom_point(aes(x=date, y=Ca, color="Ca"))+geom_line(aes(x=date, y=Ca, color="Ca"))+
  geom_point(aes(x=date, y=K, color="K orig"),data=PMF_BA_full)+geom_line(aes(x=date, y=Ca, color="Ca orig"), 
                                                                     data=PMF_BA_full)
#LD/2 para no detectados y 
#(LC-LD)/2 para lo que estan entre detectado y cuatificado.
# en los blancos  reemplazar por LD/2 