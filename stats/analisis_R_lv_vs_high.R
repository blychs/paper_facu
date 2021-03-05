library(readxl)
library(ggplot2)
library(openair)
library(tidyverse)
lv_vs_hv <- read_excel("Documents/grafico_lv_vs_hv.xlsx", 
                        sheet = "calculos intermedios")
View(lv_vs_hv) 

hist(lv_vs_hv$'HV-LV', breaks=20)
