# load libraries ####

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
# library(tidyverse)
library(patchwork)
library(readr)
library(openair)
library(lubridate)
library(readxl)
library(openairmaps)

# load dataframes ####
setwd("~/mdiaz/Documents/paper_facu/Rscripts")
Simon_mass <- read_csv("../Simon_2011_events_mass.csv",
                       col_types = cols(trace_elements = col_skip()))
Simon_mass$date <- as.POSIXct(Simon_mass$date, tz='UTC')
Simon_mass$PM2.5=rowSums(Simon_mass[,2:8])
Simon_mass$ii_perc=Simon_mass$inorganic_ions/Simon_mass$PM2.5*100
Simon_mass$om_perc=Simon_mass$organic_mass/Simon_mass$PM2.5*100
Simon_mass$gm_perc=Simon_mass$geological_minerals/Simon_mass$PM2.5*100
Simon_mass$ec_perc=Simon_mass$elemental_C/Simon_mass$PM2.5*100
Simon_mass$ss_perc=Simon_mass$salt/Simon_mass$PM2.5*100
Simon_mass$others_perc=Simon_mass$others/Simon_mass$PM2.5*100
Simon_mass$r_perc=Simon_mass$residual/Simon_mass$PM2.5*100

BA_events_testM <- read_excel("../BA_events_testMnew2.xlsx")
BA_events_testM$date <- as.POSIXct(BA_events_testM$date, tz='UTC')
BA_events_testM$Event_F <- as.factor(BA_events_testM$Event_F)

PMF_BA_full <- read_excel("../data/PMF_BA_fullv4.xlsx",
                          sheet = "CONC", col_types = c("date",
                                                        "numeric", "numeric", "numeric",
                                                        "numeric", "numeric", "numeric",
                                                        "skip", "skip", "numeric", "numeric",
                                                        "numeric", "numeric", "numeric",
                                                        "numeric", "numeric", "numeric",
                                                        "numeric", "numeric", "numeric",
                                                        "numeric", "numeric", "numeric",
                                                        "numeric", "numeric", "numeric",
                                                        "numeric", "numeric", "numeric",
                                                        "numeric", "skip", "skip", "skip",
                                                        "skip", "skip", "skip", "skip",
                                                        "skip", "skip", "skip", "skip", "skip",
                                                        "numeric", "numeric", "numeric"), na ='-999')
PMF_BA_full$date <- as.POSIXct(PMF_BA_full$date, tz='UTC')
PMF_BA_full$color <- ifelse(BA_events_testM$Event_F %in% c("SI", "SF","SO"), "BB samples", "non-BB samples")
PMF_BA_full$datefactor=factor(PMF_BA_full$date)
PMF_BA_full$date2 <- as.Date(PMF_BA_full$date, format = "%Y-%m-%d")


# prepare data ####
# Simon_mass$residual[Simon_mass$residual<0]=0
mean_data <- Simon_mass %>%
  mutate(Event_F = PMF_BA_full$color) %>%
  summarise(across(c(ii_perc, om_perc, ec_perc, gm_perc, ss_perc, others_perc, r_perc),\(x) mean(x, na.rm = TRUE))) %>%
  pivot_longer(cols = c(ii_perc, om_perc, ec_perc, gm_perc, ss_perc, others_perc, r_perc),
               names_to = "category", values_to = "value") %>%
  mutate(category = recode(category,
                           "ii_perc" = "Inorganic Ions \n 10%",
                           "om_perc" = "Organic Matter \n 57%",
                           "ec_perc" = "Elemental Carbon \n 6%",
                           "gm_perc" = "Geological Minerals \n 14%",
                           "ss_perc" = "Sea Salt \n 3%",
                           "others_perc" = "KNON \n 2%",
                           "r_perc" = "Others \n 9%"))
# Create average pie chart 
mean_data <-mean_data %>% 
  arrange(desc(category)) %>%
  mutate(percentage = value ,
         ypos = cumsum(value) - 0.5* value)

# Crear gráfico de barras de PM2.5 ####
bar_plot <- ggplot(PMF_BA_full, aes(x = date2, y = `PM2,5`, fill = color)) +
  geom_bar(stat = "identity", width = 2.5) +
  scale_fill_manual(values = c("BB samples" = "darkred", "non-BB samples" = "darkblue")) +
  scale_x_date(date_labels = "%B", date_breaks = "1 month") +
  labs(x = "", y = "PM2.5 (µg/m³)", fill = "samples") +
  theme_minimal() +theme(legend.position = "top", 
                         plot.title = element_text(size = 20, face = "bold"),  # Tamaño del título
                         text = element_text(size = 16))

bar_plot



average_pie <- ggplot(data= mean_data, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  geom_label_repel(data = mean_data,
                   aes(label = category, x = 1.3, y=ypos), size = 6, nudge_x =.9,
                   box.padding = 0.5, point.padding = 0.9,
                   segment.color = 'grey50', show.legend = TRUE,
                   direction = "y") +  # Use direction "y" to ensure labels move along the y-axis
  labs(title = "Chemical Profile") +
  theme_void() +
  theme(legend.position = "none", 
        plot.title = element_text(size = 20, face = "bold"),  # Tamaño del título
        text = element_text(size = 14)) # Hide legend
average_pie

# Combinar gráficos: barras a la izquierda, tortas a la derecha
combined_plot <- (bar_plot | (average_pie )) + 
  plot_layout(widths = c(3, 2)) 

# Mostrar gráfico combinado
print(combined_plot)
png(filename="Graphicalabstractinprogress.png", res=300, height=3*5,width = 3*13, units = "cm" )
print(combined_plot)
dev.off()

