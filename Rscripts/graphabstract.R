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

# load dataframes
setwd("~/mdiaz/Documents/paper_facu/Rscripts")
Simon_mass <- read_csv("../Simon_2011_events_mass.csv",
                       col_types = cols(trace_elements = col_skip()))
Simon_mass$date <- as.POSIXct(Simon_mass$date, tz='UTC')
Simon_mass$PM2.5=rowSums(Simon_mass[,2:8])

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

# Renaming categories in the DataFrame for pie charts
mean_data <- Simon_mass %>%
  mutate(Event_F = PMF_BA_full$color) %>%
  # group_by(Event_F) %>%
  summarise(across(c(inorganic_ions, organic_mass, elemental_C, geological_minerals, salt, others),\(x) mean(x, na.rm = TRUE))) %>%
  pivot_longer(cols = c(inorganic_ions, organic_mass, elemental_C, geological_minerals, salt, others),
               names_to = "category", values_to = "value") %>%
  mutate(category = recode(category,
                           "inorganic_ions" = "Inorganic Ions",
                           "organic_mass" = "Organic Matter",
                           "elemental_C" = "Elemental Carbon",
                           "geological_minerals" = "Geological Minerals",
                           "salt" = "Sea Salt",
                           "others" = "Others"))

# Crear gráfico de barras de PM2.5
bar_plot <- ggplot(PMF_BA_full, aes(x = date2, y = `PM2,5`, fill = color)) +
  geom_bar(stat = "identity", width = 2) +
  scale_fill_manual(values = c("BB samples" = "darkred", "non-BB samples" = "darkblue")) +
  scale_x_date(date_labels = "%B", date_breaks = "1 month") +
  labs(x = "date", y = "PM2.5 (µg/m³)", fill = "samples") +
  theme_minimal() 

bar_plot
# Create pie chart for events
mean_data <-mean_data %>%
  # filter(Event_F == "BB samples") %>%
  arrange(desc(category)) %>%
  mutate(percentage = value / sum(value) * 100,
         ypos = cumsum(value) -  value)
event_pie <- ggplot(data= mean_data, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  geom_label_repel(data = mean_data,
                   aes(label = category, x = 1.5), size = 4, nudge_x = 0.5,
                   box.padding = 0.5, point.padding = 0.5,
                   segment.color = 'grey50', show.legend = FALSE,
                   direction = "x") +  # Use direction "y" to ensure labels move along the y-axis
  labs(title = "Buenos Aires Aerosol Composition") +
  theme_void() +
  theme(legend.position = "none") # Hide legend
event_pie


# # Preparar los datos para eventos
# event_data <- mean_data %>%
#   filter(Event_F == "Event") %>%
#   arrange(desc(category)) %>%
#   mutate(percentage = value / sum(value) * 100,
#          ypos = cumsum(value) - 0.5 * value)  # Calcular posiciones para las etiquetas
# 
# # Crear gráfico de torta con geom_label_repel()
# event_pie <- ggplot(event_data) +
#   geom_bar(stat = "identity", width = 1) +
#   coord_polar("y") +
#   # Colocar las etiquetas utilizando ypos y ajustes específicos para "others" y "SS"
#   geom_label_repel(aes(y = ypos, label = category), 
#                    size = 4, 
#                    box.padding = 0.5, 
#                    point.padding = 0.5,
#                    segment.color = 'grey50', 
#                    show.legend = FALSE,
#                    direction = "y") +  # Coloca las etiquetas a lo largo del eje Y
#   labs(title = "BB samples") +
#   theme_void() +
#   theme(legend.position = "none")  # Ocultar leyenda
# 
# event_pie
# 
# 
# 
# 
# # Crear gráfico de pastel para eventos sin etiquetas
# no_event_pie <- ggplot(mean_data %>% filter(Event_F == "No Event"), aes(x = "", y = value, fill = category)) +
#   geom_bar(stat = "identity", width = 1) +
#   coord_polar("y") +
#   # Etiquetas dentro del gráfico para categorías que caben
#   geom_text(data = filter(mean_data %>% filter(Event_F == "No Event")),
#             aes(label = category), 
#             position = position_stack(vjust = 0.5), size = 4, 
#             hjust = 0.5) +
#   labs(title = "No Event") +
#   theme_void() +
#   theme(legend.position = "none") # Ocultar leyenda
# 
# # Combine plots: bars on the left, pies on the right
# combined_plot <- (bar_plot | (event_pie / no_event_pie))
# combined_plot <- (bar_plot | (event_pie ))
# # Show combined plot
# print(combined_plot)


# Combinar gráficos: barras a la izquierda, tortas a la derecha
combined_plot <- (bar_plot | (event_pie ))

# Mostrar gráfico combinado
print(combined_plot)

