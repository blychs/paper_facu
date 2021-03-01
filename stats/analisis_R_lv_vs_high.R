library("ggplot2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("openair", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")

lv_vs_hv <- read_excel("Documents/graficos_lv_vs_hv.xlsx", 
                        sheet = "calculos intermedios")
View(graficos_lv_vs_hv) 

hist(graficos_lv_vs_hv$'HV-LV', breaks=20)
