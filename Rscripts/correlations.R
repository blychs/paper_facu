# 00 Load libraries ####
library(openair)
library(ggplot2)
library(lubridate)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(openairmaps)
pathgraphs="Figures"
MWNH4=18.04
MWNO3=62.0049
MWSO4=96.06
MWK=39.0898
MWCa=40.08
MWMg=24.31
MWNa=23
MWCl=35.45

setwd("~/mdiaz/Documents/paper_facu/Rscripts")
BA_events_testM <- read_excel("../data/BA_events_testMnew2.xlsx")
BA_events_testM$date <- as.POSIXct(BA_events_testM$date, tz='UTC')
BA_events_testM$Event_F <- as.factor(BA_events_testM$Event_F)
BA_events_testM$Event_F=factor(BA_events_testM$Event_F)

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
                                                        "numeric", "numeric", "numeric"),
                          na ='-999')

PMF_BA_full=dplyr::rename(PMF_BA_full, EC = 'C Elemental',TC = 'C Total',OC = 'C Orgánico')
# PMF_BA_full_so=PMF_BA_full[PMF_BA_full$date<=as.POSIXct('2019-05-23') | PMF_BA_full$date>=as.POSIXct('2019-06-02'),]
PMF_BA_full$OC_EC = PMF_BA_full$OC/PMF_BA_full$EC
PMF_BA_full$OC_K=PMF_BA_full$OC/PMF_BA_full$K
PMF_BA_full$KOC=PMF_BA_full$K/PMF_BA_full$OC

# 01 Matriz de correlaciones ####
png("../images/corrcoefmaps.png", res = 300, height = 3000, width = 3000)
corPlot(PMF_BA_full, dendrogram = TRUE,method = "pearson", main= "R pearson")
dev.off()
png("../images/corrcoefmaps_S.png", res = 300, height = 3000, width = 3000)
corPlot(PMF_BA_full, dendrogram = TRUE,method = "spearman", main ="R spearman, matriz completa")
dev.off()

colMeans(PMF_BA_full[,sapply(PMF_BA_full, is.numeric)], na.rm = TRUE)


#GA ####
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

# Renaming categories in the DataFrame
mean_data <- Simon_mass %>%
  mutate(Event_F = PMF_BA_full$color) %>%
  group_by(Event_F) %>%
  summarise(across(c(inorganic_ions, organic_mass, elemental_C, geological_minerals, salt, others), mean, na.rm = TRUE)) %>%
  pivot_longer(cols = c(inorganic_ions, organic_mass, elemental_C, geological_minerals, salt, others),
               names_to = "category", values_to = "value") %>%
  mutate(category = recode(category,
                           "inorganic_ions" = "II",
                           "organic_mass" = "OM",
                           "elemental_C" = "EC",
                           "geological_minerals" = "GM",
                           "salt" = "SS",
                           "others" = "others"))

# Define categories to be shown outside the pie chart

outside_labels <- c("EC","SS", "others")


library(ggplot2)
library(ggrepel)
library(tidyverse)
library(patchwork)
# Crear gráfico de barras de PM2.5
PMF_BA_full$color <- ifelse(PMF_BA_full$Event_F %in% c("SI", "SF","SO"), "BB samples", "non-BB samples")

PMF_BA_full$date=factor(PMF_BA_full$date)
bar_plot <- ggplot(PMF_BA_full, aes(x = date, y = `PM2,5`, fill = color)) +
  geom_bar(stat = "identity", width = 0.95) +
  scale_fill_manual(values = c("BB samples" = "red", "non-BB samples" = "blue")) +
  scale_x_discrete(breaks = levels(PMF_BA_full$date)[seq(1, length(levels(PMF_BA_full$date)), by = 4)]) +
  labs(x = "Date", y = "PM2.5 (µg/m³)", fill = "Samples") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# # Create pie chart for events
# event_data <-mean_data %>%
#   filter(Event_F == "Event") %>%
#   arrange(desc(category)) %>%
#   mutate(percentage = value / sum(value) * 100,
#          ypos = cumsum(value) - 0.5 * value) 
# event_pie <- ggplot(event_data, aes(x = "", y = value, fill = category)) +
#   geom_bar(stat = "identity", width = 1) +
#   coord_polar("y") +
#   geom_label_repel(data = event_data,
#                    aes(label = category, x = 1.1), size = 4, nudge_x = 0.5,
#                    box.padding = 0.5, point.padding = 0.5,
#                    segment.color = 'grey50', show.legend = FALSE,
#                    direction = "y") +  # Use direction "y" to ensure labels move along the y-axis
#   labs(title = "BB samples") +
#   theme_void() +
#   theme(legend.position = "none") # Hide legend
# event_pie
library(ggplot2)
library(dplyr)
library(ggrepel)

# Preparar los datos para eventos
event_data <- mean_data %>%
  filter(Event_F == "Event") %>%
  arrange(desc(category)) %>%
  mutate(percentage = value / sum(value) * 100,
         ypos = cumsum(value) - 0.5 * value)  # Calcular posiciones para las etiquetas

# Crear gráfico de torta con geom_label_repel()
event_pie <- ggplot(event_data, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  
  # Colocar las etiquetas utilizando ypos y ajustes específicos para "others" y "SS"
  geom_label_repel(aes(y = ypos, label = category), 
                   size = 4, 
                   nudge_x = ifelse(event_data$category %in% c("others", "SS"), 1, 0.5),  # Ajuste mayor para others y SS
                   nudge_y = ifelse(event_data$category %in% c("others", "SS"), 0.2, 0),  # Ajustar verticalmente si es necesario
                   box.padding = 0.5, 
                   point.padding = 0.5,
                   segment.color = 'grey50', 
                   show.legend = FALSE,
                   direction = "y") +  # Coloca las etiquetas a lo largo del eje Y
  
  labs(title = "BB samples") +
  theme_void() +
  theme(legend.position = "none")  # Ocultar leyenda

event_pie




# Crear gráfico de pastel para eventos sin etiquetas
no_event_pie <- ggplot(mean_data %>% filter(Event_F == "No Event"), aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  # Etiquetas dentro del gráfico para categorías que caben
  geom_text(data = filter(mean_data %>% filter(Event_F == "No Event")),
            aes(label = category), 
            position = position_stack(vjust = 0.5), size = 4, 
            hjust = 0.5) +
 labs(title = "No Event") +
  theme_void() +
  theme(legend.position = "none") # Ocultar leyenda

# Combine plots: bars on the left, pies on the right
combined_plot <- (bar_plot | (event_pie / no_event_pie))

# Show combined plot
print(combined_plot)


# Combinar gráficos: barras a la izquierda, tortas a la derecha
combined_plot <- (bar_plot / (event_pie | no_event_pie))

# Mostrar gráfico combinado
print(combined_plot)


# 03 Mergeo con eventos ####
PMF_BA_full=merge(PMF_BA_full,BA_events_testM, by="date")
PMF_BA_full$NH4molar=PMF_BA_full$NH4/MWNH4
PMF_BA_full$NO3molar=PMF_BA_full$NO3/MWNO3
PMF_BA_full$SO4molar=PMF_BA_full$SO4/MWSO4
PMF_BA_full$Namolar=PMF_BA_full$`Na sol`/MWNa
PMF_BA_full$Kmolar=PMF_BA_full$K/MWK
PMF_BA_full$NH4Simon=PMF_BA_full$NH4molar
PMF_BA_full$NO3SO4Simon=PMF_BA_full$NO3molar+2*PMF_BA_full$SO4molar
PMF_BA_full$NeutralizedSimon=PMF_BA_full$NH4molar-(PMF_BA_full$NO3molar+2*PMF_BA_full$SO4molar)
PMF_BA_full$NH4NO3 = (PMF_BA_full$NH4molar)/(PMF_BA_full$NO3molar)
PMF_BA_full$NH4SO4 = (PMF_BA_full$NH4molar)/(PMF_BA_full$SO4molar)
PMF_BA_full$NH42SO4NO3 = (PMF_BA_full$NH4molar)/((2*PMF_BA_full$SO4molar)+(PMF_BA_full$NO3molar))
PMF_BA_full$NH4SO4NO3 = (PMF_BA_full$NH4molar)/((PMF_BA_full$SO4molar)+(PMF_BA_full$NO3molar))
PMF_BA_full$KSO4 = (PMF_BA_full$Kmolar)/((PMF_BA_full$SO4molar))
PMF_BA_full$NO3SO4 = (PMF_BA_full$NO3molar)/((PMF_BA_full$SO4molar))
PMF_BA_full$ClNa = (PMF_BA_full$Cl)/((PMF_BA_full$`Na sol`))
PMF_BA_full$ClNamolar = (PMF_BA_full$Cl/MWCl)/((PMF_BA_full$`Na sol`/MWNa))
PMF_BA_full$Na_perc=PMF_BA_full$`Na sol`/PMF_BA_full$`PM2,5`
# conteo con eventos ####
event_levels <- c("SI", "SF")
summary(PMF_BA_full[PMF_BA_full$Event_F %in% event_levels,])
summary(PMF_BA_full[PMF_BA_full$Event_F %in% c("no"),])

PMFBAsorted<-PMF_BA_full[order(PMF_BA_full$`PM2,5`),]
PMFBAsorted<-PMFBAsorted[,c("date", "PM2,5", "Event_F")]
# Count the number of events 
event_levels <- c("SI", "SF")
event_MAM <- sum(table(BA_events_testM$Event_F[month(BA_events_testM$date)>=3&month(BA_events_testM$date)<=5])[event_levels])
event_JJA <- sum(table(BA_events_testM$Event_F[month(BA_events_testM$date)>=6&month(BA_events_testM$date)<=8])[event_levels])
event_SON <- sum(table(BA_events_testM$Event_F[month(BA_events_testM$date)>=9&month(BA_events_testM$date)<=11])[event_levels])
event_DJF <- sum(table(BA_events_testM$Event_F[month(BA_events_testM$date)==12|month(BA_events_testM$date)<=2])[event_levels])
print(paste("Cant. de eventos: MAM", event_MAM, ", JJA:", event_JJA, ", SON:", event_SON, ", DJF:", event_DJF))

event_JAS <- sum(table(BA_events_testM$Event_F[month(BA_events_testM$date)>=7&month(BA_events_testM$date)<=9])[event_levels])
print(event_JAS)


ggplot(PMF_BA_full)+geom_line(aes(x=date,y=ClNamolar))+geom_abline(slope=0, intercept = 1.16,linetype="dashed")
#  04 Iones ####
ggplot(PMF_BA_full)+
  geom_line(aes(x=date, y=Kmolar,color="K"))+
  geom_line(aes(x=date, y=NH4molar,color="NH4"))+
  geom_line(aes(x=date, y=NO3molar,color="NO3"))+
  geom_line(aes(x=date, y=SO4molar,color="SO4"))+
  geom_point(data=PMF_BA_full[PMF_BA_full$Event_F %in% c("SI","SF"),],aes(x=date, y=NH4molar),color="black",show.legend = FALSE)

ggplot(PMF_BA_full)+
  geom_point(aes(x=NO3molar, y=(NH4molar/SO4molar-2)*SO4molar,color="EXCEDENCIA"))
i1=1
i2=2
library(ggpubr)
ratios <- ggplot(PMF_BA_full) + 
  geom_line(aes(x=date, y=NH4SO4, color="NH4/SO4[molar]"))+
  geom_line(aes(x=date, y=NH42SO4NO3, color="NH4/(2SO4+NO3)[molar]"))+
  geom_abline(aes(slope=0,intercept=i1), linetype="dashed")+
  geom_abline(aes(slope=0,intercept=i2), color='gray20' ,linetype="dashed")+
  labs(y="molar ratio",x="")+  
  geom_point(data = PMF_BA_full[PMF_BA_full$Event_F %in% c("SI"),], 
             aes(x=date,y=NH4SO4, fill=OrigenTag), shape=21, size=2)+
  scale_color_discrete(name = "Molar Ratios", na.translate = F)+
  scale_fill_discrete(name = "Origin", na.translate = F)+
  guides(color = guide_legend(nrow = 2),fill = guide_legend(nrow = 3)) + 
  scale_y_continuous(
    breaks = function(x) unique(c(i1,i2, scales::pretty_breaks()(x))),
    labels = function(x) format(x, nsmall = 0, digits = 2)
  )+theme_bw()+theme(legend.position="top")

#NH4 NO3 SO4 y eventos
ions <- ggplot(PMF_BA_full) + 
    geom_line(aes(x=date, y=NH4, color="NH4"))+
    geom_line(aes(x=date, y=NO3, color="NO3"))+
    geom_line(aes(x=date, y=SO4, color="SO4"))+
    labs(y=expression(paste("concentrations [", mu, "g/m"^3, "]")),x="")+  
    geom_point(data = PMF_BA_full[PMF_BA_full$Event_F %in% c("SI"),], 
               aes(x=date,y=NO3, fill=OrigenTag), shape=21, size=2)+
    scale_color_discrete(name = "Ions", na.translate = F)+
    scale_fill_discrete(name = "Origin", na.translate = F)+
    guides(color = guide_legend(nrow = 2),fill = guide_legend(nrow = 3)) + 
    scale_y_continuous(
      breaks = function(x) unique(c(i1,i2, scales::pretty_breaks()(x))),
      labels = function(x) format(x, nsmall = 0, digits = 2)
    )+theme_bw()+theme(legend.position="top")

combined_plot <- ggarrange(ratios, ions, ncol = 1, nrow = 2, common.legend = TRUE, legend = "top")

# Mostrar el gráfico combinado
print(combined_plot)

# grafico 8 ####
# Crear los gráficos por separado
library(patchwork)
plot1 <- ggplot(PMF_BA_full) + 
  geom_line(aes(x = date, y = NH4SO4, color = "NH4/SO4"), linewidth = 1.1 ) +
  geom_line(aes(x = date, y = NH42SO4NO3, color = "NH4/(2SO4+NO3)"), linewidth = 1.1 ) +
  geom_abline(aes(slope = 0, intercept = i1), linetype = "dashed") +
  geom_abline(aes(slope = 0, intercept = i2), color = 'gray20', linetype = "dashed") +
  labs(y = "molar ratio", x = "") +  
  # geom_point(data = PMF_BA_full[PMF_BA_full$Event_F %in% c("SI","SF","SO"), ], 
  #            aes(x = date, y = NH4SO4, fill = OrigenTag), shape = 21, size = 2) +
  geom_point(data = PMF_BA_full[PMF_BA_full$Event_F %in% c("SI","SF","SO"), ], 
             aes(x = date, y = NH4SO4), fill = "black", shape = 21, size = 2) +
  scale_color_discrete(name = "molar ratios", na.translate = FALSE) +
  # scale_fill_discrete(name = "SE-origin", na.translate = FALSE) +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) + 
  scale_y_continuous(
    breaks = function(x) unique(c(i1, i2, scales::pretty_breaks()(x))),
    labels = function(x) format(x, nsmall = 0, digits = 2)
  ) + 
  theme_bw() + 
  theme(
    legend.position = "top",
    axis.title = element_text(size = 14),  # Tamaño del título del eje
    axis.text = element_text(size = 14),   # Tamaño del texto de los ticks del eje
    legend.text = element_text(size = 14), # Tamaño del texto de la leyenda
    legend.title = element_text(size = 14), # Tamaño del título de la leyenda
    plot.margin=unit(c(0,0,0,0), "pt")
  ) +
  annotate("text", x = min(PMF_BA_full$date), y = Inf, label = "(a)", 
           vjust = 1.5, hjust = 1.4, size = 5, fontface = "bold")+
  scale_x_datetime(expand = c(0, 0),
                   date_labels ="%m-%Y")

plot2 <- ggplot(PMF_BA_full) + 
  geom_line(aes(x = date, y = NH4, color = "NH4"), linewidth = 1.1 ) +
  geom_line(aes(x = date, y = NO3, color = "NO3"), linewidth = 1.1 ) +
  geom_line(aes(x = date, y = SO4, color = "SO4"), linewidth = 1.1 ) +
  labs(y = expression(paste("concentration [", mu, "g/m"^3, "]")), x = "") +  
  # geom_point(data = PMF_BA_full[PMF_BA_full$Event_F %in% c("SI","SF","SO"), ], 
  #            aes(x = date, y = NO3, fill = OrigenTag), shape = 21, size = 2) +
  geom_point(data = PMF_BA_full[PMF_BA_full$Event_F %in% c("SI","SF","SO"), ], 
             aes(x = date, y = NO3),fill="black", shape = 21, size = 2) +
  scale_color_discrete(name = "ions concentration", na.translate = FALSE) +
  # scale_fill_discrete(name = "SE-origin", na.translate = FALSE) +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 2)) + 
  scale_y_continuous(
    breaks = function(x) unique(c(i1, i2, scales::pretty_breaks()(x))),
    labels = function(x) format(x, nsmall = 0, digits = 2)
  ) + 
  theme_bw() + 
  theme(
    legend.position = "top",
    axis.title = element_text(size = 14),  # Tamaño del título del eje
    axis.text = element_text(size = 14),   # Tamaño del texto de los ticks del eje
    legend.text = element_text(size = 14), # Tamaño del texto de la leyenda
    legend.title = element_text(size = 14), # Tamaño del título de la leyenda
    plot.margin=unit(c(0,0,0,0), "pt")
  ) +
  annotate("text", x = min(PMF_BA_full$date), y = Inf, label = "(b)", 
           vjust = 1.5, hjust = 1.4, size = 5, fontface = "bold")+
  scale_x_datetime(expand = c(0, 0),
                   date_labels ="%m-%Y")

# Usar patchwork para combinar los gráficos y conservar las tres leyendas
combined_plot <- plot1 / plot2 + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = "top",  legend.box.just = "top",legend.box = "vertical") &
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 2))

# Mostrar el gráfico combinado
print(combined_plot)
# hasta aca ####

png(paste0("../",pathgraphs,"/ionesratiosblack.png"), width = 550 * 5, height = 700* 5, res = 300)
print(combined_plot)
dev.off()

# molar ratio tesis ####
library(patchwork)
plot1 <- ggplot(PMF_BA_full) + 
  geom_line(aes(x = date, y = NH4SO4, color = "NH4/SO4"), linewidth = 1.1 ) +
  geom_line(aes(x = date, y = NH42SO4NO3, color = "NH4/(2SO4+NO3)"), linewidth = 1.1 ) +
  geom_abline(aes(slope = 0, intercept = i1), linetype = "dashed") +
  geom_abline(aes(slope = 0, intercept = i2), color = 'gray20', linetype = "dashed") +
  labs(y = "relación en moles", x = "") +  
  geom_point(data = PMF_BA_full[PMF_BA_full$Event_F %in% c("SI","SF","SO"), ], 
             aes(x = date, y = NH4SO4), fill = "black", shape = 21, size = 2) +
  scale_color_manual(name = "relación en moles", 
                     values = c("NH4/SO4" = "blue", "NH4/(2SO4+NO3)" = "red"),
                     labels = c(expression(NH[4]^"+" / SO[4]^"2-"), 
                                expression(NH[4]^"+" / (2*SO[4]^"2-" + NO[3]^"-")))) +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) + 
  scale_y_continuous(
    breaks = function(x) unique(c(i1, i2, scales::pretty_breaks()(x))),
    labels = function(x) format(x, nsmall = 0, digits = 2)
  ) + 
  theme_bw() + 
  theme(
    legend.position = "top",
    axis.title = element_text(size = 14),  # Tamaño del título del eje
    axis.text = element_text(size = 14),   # Tamaño del texto de los ticks del eje
    legend.text = element_text(size = 14), # Tamaño del texto de la leyenda
    legend.title = element_text(size = 14), # Tamaño del título de la leyenda
    plot.margin=unit(c(0,0,0,0), "pt")
  ) +
  annotate("text", x = min(PMF_BA_full$date), y = Inf, label = "(a)", 
           vjust = 1.5, hjust = 1.4, size = 5, fontface = "bold")+
  scale_x_datetime(expand = c(0, 0),
                   date_labels ="%m-%Y")

plot2 <- ggplot(PMF_BA_full) + 
  geom_line(aes(x = date, y = NH4, color = "NH4"), linewidth = 1.1 ) +
  geom_line(aes(x = date, y = NO3, color = "NO3"), linewidth = 1.1 ) +
  geom_line(aes(x = date, y = SO4, color = "SO4"), linewidth = 1.1 ) +
  labs(y = expression(paste("concentración [", mu, "g/m"^3, "]")), x = "") +  
  geom_point(data = PMF_BA_full[PMF_BA_full$Event_F %in% c("SI","SF","SO"), ], 
             aes(x = date, y = NO3), fill="black", shape = 21, size = 2) +
  scale_color_discrete(name = "iones", 
                     labels = c(expression(NH[4]^"+"), expression(NO[3]^"-"), expression(SO[4]^"2-"))) +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 2)) + 
  scale_y_continuous(
    breaks = function(x) unique(c(i1, i2, scales::pretty_breaks()(x))),
    labels = function(x) format(x, nsmall = 0, digits = 2)
  ) + 
  theme_bw() + 
  theme(
    legend.position = "top",
    axis.title = element_text(size = 14),  # Tamaño del título del eje
    axis.text = element_text(size = 14),   # Tamaño del texto de los ticks del eje
    legend.text = element_text(size = 14), # Tamaño del texto de la leyenda
    legend.title = element_text(size = 14), # Tamaño del título de la leyenda
    plot.margin=unit(c(0,0,0,0), "pt")
  ) +
  annotate("text", x = min(PMF_BA_full$date), y = Inf, label = "(b)", 
           vjust = 1.5, hjust = 1.4, size = 5, fontface = "bold")+
  scale_x_datetime(expand = c(0, 0),
                   date_labels ="%m-%Y")

# Usar patchwork para combinar los gráficos y conservar las tres leyendas
combined_plot <- plot1 / plot2 + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = "top",  legend.box.just = "top",legend.box = "vertical") &
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 2))

# Mostrar el gráfico combinado
print(combined_plot)
# hasta aca ####

png(paste0("../",pathgraphs,"/ionesratios_tesis.png"), width = 550 * 5, height = 700* 5, res = 300)
print(combined_plot)
dev.off()
# neutralized ####

png("../images/neutralizedSimon.png", res = 300, height = 3000, width = 3000)
ggplot(PMF_BA_full)+geom_line(aes(x=date,y=NeutralizedSimon),color="magenta")+geom_abline(aes(slope=0,intercept=0),linetype="dashed")
dev.off()

keys = names(PMF_BA_full)[!names(PMF_BA_full) %in% c("Event", "date", "date mal", "ID" , "Event_M","Razon","Event_F","Local","Tag","Origen" ,"OrigenTag",  "Observaciones" , "Lote")]
print(keys)
colors = setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(keys)), keys)
pathgraphs = "Figures"
while (!is.null(dev.list())) dev.off()

for (title in keys) {
  tryCatch({
    safe_title = make.names(title)
    output_path = file.path(getwd(), "../Figures", paste0("serietemporal_", safe_title, ".png"))
    p <- ggplot(PMF_BA_full) +
      geom_line(aes(x = date, y = !!sym(title)), color = colors[title]) +
      ggtitle(paste("Plot for", title))  # Añadir título a la gráfica para depuración
    png(output_path, width = 590 * 3, height = 592 * 3, res = 300)
    print(p)
    dev.off()
  }, error = function(e) {
    error_message <- paste("Error en el gráfico de", title, ": ", e$message)
    message(error_message)
  })
}

# 04 b Periodos Iones ####

summary(selectByDate(PMF_BA_full,end = "2019-05-25"))
summary(selectByDate(PMF_BA_full,start = "2019-05-25",end = "2019-10-01"))
summary(selectByDate(PMF_BA_full,start = "2019-10-01"))

#Hasta Mayo
# 
pruebacut=selectByDate(PMF_BA_full, end="2019-05-25")
ggplot(pruebacut)+
  geom_point(aes(x=NH4molar, y=NO3molar))+
  geom_abline(aes(slope=1, intercept=0),linetype="dashed")
ggplot(pruebacut)+
  geom_point(aes(x=NH4molar, y=SO4molar))+
  geom_abline(aes(slope=1, intercept=0),linetype="dashed")
ggplot(pruebacut)+
  geom_point(aes(x=NH4molar, y=NO3molar))+
  geom_abline(aes(slope=1, intercept=0),linetype="dashed")

corPlot(PMF_BA_full[,c("Cl","SO4","NH4","NO3","K","Na sol","Mg")], dendrogram = TRUE, 
        method = "pearson")
png("../images/corplothastaMayo.png", res = 300, height = 3000, width = 3000)
corPlot(pruebacut[,c("SO4","NH4","NO3","K","Na sol","Mg")], dendrogram = TRUE, 
        method = "pearson", main ="Hasta Mayo")
dev.off()

ggplot(pruebacut) + geom_point(aes(x=Kmolar,y=NO3molar))+geom_point(data=pruebacut[pruebacut$Event_F %in% c("SI","SF"),],aes(x=Kmolar, y=NO3molar),color="red",show.legend = FALSE)
#Octubre en adelante
pruebacut=selectByDate(PMF_BA_full, start="2019-10-01")

ggplot(PMF_BA_full,aes(x=NO3molar,y=NH4molar))+geom_point()+xlim(0,0.18)+ylim(0,0.18)+
  stat_summary(fun.data= mean_cl_normal) + 
  geom_smooth(method='lm',se = F)
  coefs <- coef(lm(NH4molar~NO3molar, data = PMF_BA_full))
coefs

ggplot(PMF_BA_full,aes(x=2*SO4molar,y=NH4molar))+geom_point()+
  stat_summary(fun.data= mean_cl_normal) + 
  geom_smooth(method='lm',se = F)
PMF_BA_full$SO4molar_doble <- 2 * PMF_BA_full$SO4molar

# Ajustar el modelo de regresión
modelo <- lm(NH4molar ~ SO4molar_doble, data = PMF_BA_full, na.action = na.omit)
coefs <- coef(lm(NH4molar~SO4molar_doble, data = PMF_BA_full, na.action = na.omit))
coefs

ggplot(pruebacut) + geom_point(aes(x=Kmolar,y=NO3molar))+geom_point(data=pruebacut[pruebacut$Event_F %in% c("SI","SF"),],aes(x=Kmolar, y=NO3molar),color="red",show.legend = FALSE)

ggplot(pruebacut)+
  geom_point(aes(x=NH4molar, y=NO3molar))+xlim(0,0.03)+ylim(0,0.03)+
  geom_abline(aes(slope=1, intercept=0),linetype="dashed")

ggplot(pruebacut)+
  geom_point(aes(y=Namolar, x=SO4molar))+
  geom_abline(aes(slope=1, intercept=0),linetype="dashed")+xlim(0,0.03)+ylim(0,0.03)

ggplot(pruebacut)+
  geom_point(aes(x=NH4molar, y=SO4molar))+
  geom_abline(aes(slope=1, intercept=0),linetype="dashed")+xlim(0,0.03)+ylim(0,0.03)

ggplot(pruebacut)+
  geom_point(aes(x=Kmolar, y=NO3molar))+
  geom_abline(aes(slope=1, intercept=0),linetype="dashed")+xlim(0,0.03)+ylim(0,0.03)

png("../images/corplotdesdeOctubre.png", res = 300, height = 3000, width = 3000)
corPlot(pruebacut[,c("SO4","NH4","NO3","K","Na sol","Mg")], dendrogram = TRUE, 
        method = "pearson", main ="desdeoctubre")
dev.off()
png("../images/NaSO4desdeOctubre.png", res = 300, height = 3000, width = 3000)
ggplot(pruebacut) + geom_point(aes(x=`Na sol` ,y=SO4))
dev.off()

# mayo a octubre
pruebacut=selectByDate(PMF_BA_full, start="2019-05-25",end="2019-10-01")
ggplot(pruebacut)+
  geom_point(aes(x=NH4molar, y=NO3molar))+
  geom_abline(aes(slope=1, intercept=0),linetype="dashed")
ggplot(pruebacut)+
  geom_point(aes(x=NH4molar, y=SO4molar))+
  geom_abline(aes(slope=0.5, intercept=0),linetype="dashed")

png("../images/corplotMayoaOctubre.png", res = 300, height = 3000, width = 3000)
corPlot(pruebacut[,c("SO4","NH4","NO3","K","Na sol","Mg")], dendrogram = TRUE, 
        method = "pearson", main ="Mayo a octubre")
dev.off()
ggplot(pruebacut) + geom_point(aes(x=`Na sol` ,y=SO4))
ggplot(pruebacut) + geom_point(aes(x=Kmolar,y=NO3molar))+geom_point(data=pruebacut[pruebacut$Event_F %in% c("SI","SF"),],aes(x=Kmolar, y=NO3molar),color="red",show.legend = FALSE)

# 04 c Ratios IONES ####


keys =  c("NH4NO3", "NH4SO4", "NH4SO4NO3", "NH42SO4NO3" , "KSO4","NO3SO4")
print(keys)
colors = setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(keys)), keys)
pathgraphs = "Figures"
while (!is.null(dev.list())) dev.off()

for (title in keys) {
  tryCatch({
    safe_title = make.names(title)
    output_path = file.path(getwd(), "../Figures", paste0("serietemporal_ratio_", safe_title, ".png"))
    p <- ggplot()+geom_line(data=PMF_BA_full, aes(x=date,y=!!sym(title)))+
      geom_abline(slope=0,intercept=1, linetype=3, color="black")+
      geom_abline(slope=0,intercept=1.5, linetype=3, color="magenta")+
      geom_abline(slope=0,intercept=2, linetype=3, color="blue")+
      geom_point(data=PMF_BA_full[PMF_BA_full$Event_F %in% c("SI","SF"),],aes(x=date, y=!!sym(title)),color="red",show.legend = FALSE)
      ggtitle(paste("Plot for ratio", title))  # Añadir título a la gráfica para depuración
    png(output_path, width = 590 * 3, height = 592 * 3, res = 300)
    print(p)
    dev.off()
  }, error = function(e) {
    error_message <- paste("Error en el gráfico de", title, ": ", e$message)
    message(error_message)
  })
}


ggplot(PMF_BA_full)+
  geom_point(aes(x=NH4molar, y=NO3molar))+
  geom_abline(aes(slope=1, intercept=0),linetype="dashed")

# 05 Metales ####
corPlot(PMF_BA_full[,c("Cu","Mo","Sb","Zn","Ni","Pb")], dendrogram = TRUE,method = "pearson", main ="R pearson, matriz completa")
ggplot(PMF_BA_full) + geom_point(aes(x=Pb ,y=Sb))
ggplot(PMF_BA_full) + geom_point(aes(x=Cu ,y=Sb))

png(file = "~/mdiaz/Documents/paper_facu/images/MoSb_scatterplot.png", res=300, height= 300*4, width = 4*400)
ggplot(PMF_BA_full) + geom_point(aes(x=Mo ,y=Sb))+ xlim(0,0.003) 
dev.off()

ggplot(PMF_BA_full) + geom_point(aes(x=Zn ,y=Sb))
ggplot(PMF_BA_full) + geom_point(aes(x=Cu ,y=Mo)) + ylim(0,0.003)


scatterPlot(PMF_BA_full, x="Ni", y="V")
ggplot(PMF_BA_full)+ geom_point(aes(x=Ni, y=V))+geom_abline(slope = 0.7)+geom_abline(slope = 0.15)+geom_abline(slope =2.8)

# 06 geological minerals ####
corPlot(PMF_BA_full[,c("Al","Ba", "Ca","Fe", "Mn", "Mg", "Ti","Sb")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson, sin outliers")
corPlot(PMF_BA_full[,c( "Ca","Ba",  "Mg", "Ti","Sb")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson, sin outliers")
corPlot(PMF_BA_full[,c( "Fe","Al",  "Mn", "Ti")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson, sin outliers")

corPlot(PMF_BA_full[,c( "Fe","Al",  "Mn", "Ti")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson, sin outliers")
corPlot(PMF_BA_full[,c( "Fe","Al","Ba","Cu","Pb","Sb","Zn", "As")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson, traffic related")
corPlot(PMF_BA_full[,c("Zn", "Cu","Pb")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson, traffic related")


# 07 nss ####

PMF_BA_full$nssK=PMF_BA_full$K-0.6*PMF_BA_full$Fe-0.037*PMF_BA_full$`Na sol`
PMF_BA_full$nssK_OC=PMF_BA_full$nssK/PMF_BA_full$OC
PMF_BA_full$nssK_EC=PMF_BA_full$nssK/PMF_BA_full$EC
PMF_BA_full$OC_EC=PMF_BA_full$OC/PMF_BA_full$EC
corPlot(PMF_BA_full[,c( "Na sol","Cl",  "nssK", "OC","EC")], 
        dendrogram = TRUE, main= "R salt")



corPlot(PMF_BA_full[,c( "Na sol","Cl",  "NH4", "NO3","SO4")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson")


# 08 biomass burning ####
corPlot(PMF_BA_full[,c( "OC","EC","K","NO3","NH4","SO4")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson")
corPlot(PMF_BA_full[PMF_BA_full$Event_F %in% c("SI","SF","SN"),c( "OC","EC","K","NO3","NH4","SO4")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson")
corPlot(PMF_BA_full[PMF_BA_full$Event_F %in% c("no"),c( "OC","EC","K","NO3","NH4","SO4")], 
        dendrogram = TRUE,method = "pearson", main= "R pearson")


corPlot(selectByDate(PMF_BA_full, start = "2019-07-25", end = "2019-09-25"),
        pollutants = c("NO3","NH4","K","EC","OC","SO4"),dendrogram = TRUE)
corPlot(PMF_BA_full,dendrogram = TRUE)

PMF_BA_full$neutralization=PMF_BA_full$NH4/(PMF_BA_full$NO3+PMF_BA_full$SO4)
ggplot(PMF_BA_full)+
  geom_line(aes(x=date, y=neutralization))+
  geom_line(aes(x=date, y=K,color="K"))+
  geom_line(aes(x=date, y=SO4,color="Cu"))



# 09 Time Variations ####
TVmean=timeVariation(PMF_BA_full,pollutant = "PM2,5")
print(TVmean$plot$month)

timePlot(PMF_BA_full,pollutant = "PM2,5")
timePlot(PMF_BA_full,pollutant = "PM2,5", avg.time = "month")

timePlot(PMF_BA_full,pollutant = "KOC")
ggplot(PMF_BA_full)+
  geom_line(aes(x=date, y=KOC,color="K/OC"))+
  geom_point(data=PMF_BA_full[PMF_BA_full$Event_F %in% c("SI","SF"),],aes(x=date, y=KOC),color="black",show.legend = FALSE)

ggplot(PMF_BA_full)+
  geom_line(aes(x=date, y=KOC,color="K/OC"))+
  geom_point(data=PMF_BA_full[PMF_BA_full$Event_F %in% c("SI","SF","SN","SP"),],aes(x=date, y=KOC,color=OrigenTag),show.legend = FALSE)+
  geom_point(data=PMF_BA_full[PMF_BA_full$Event_F %in% c("no"),],aes(x=date, y=KOC,color=Event_F),show.legend = FALSE)

X2C5f_base <- read_excel("~/Documents/paper_facu/data/2C5f_base.xlsx", sheet = "Contributions", range = "B108:H208")

ggplot(PMF_BA_full)+
  geom_line(aes(x=date, y=KOC,color="K/OC")) + 
  geom_line(aes(x=date, y=`PM2,5`/50,color="PM2.5")) +
  geom_line(data=X2C5f_base, aes(x=PMF_BA_full$date,y=`Factor 5`/10))+
  geom_point(data=PMF_BA_full[PMF_BA_full$Event_F %in% c("SI","SF","SN","SP"),],aes(x=date, y=`PM2,5`/50,color=OrigenTag),show.legend = FALSE)
  # geom_point(data=PMF_BA_full[PMF_BA_full$Event_F %in% c("no"),],aes(x=date, y=KOC,color=Event_F),show.legend = FALSE)

print(PMF_BA_full$date[PMF_BA_full$KOC>0.2])

# gases ####
gases20192020day=timeAverage(gases20192020, avg.time = "day")

gases20192020$date=gases20192020$day+gases20192020$Hora
ggplot(selectByDate(PMF_BA_full,start = "2019-05-15",end = "2019-06-15"))+
  geom_line(aes(x=date,y=SO4,color="SO4"))+
  geom_line(aes(x=date,y=NO3,color="NO3"))+
  geom_line(aes(x=date,y=`Na sol`,color="Na"))+
  geom_line(aes(x=date,y=Cl,color="Cl"))+
  geom_line(aes(x=date,y=NH4,color="NH4"))

ggplot(PMF_BA_full)+geom_line(aes(x=date,y=SO4,color="SO4"))+geom_point(aes(x=date,y=SO4,color=OrigenTag))
ggplot(PMF_BA_full)+geom_line(aes(x=date,y=NO3,color="NO3"))+geom_point(aes(x=date,y=NO3,color=OrigenTag))
ggplot(PMF_BA_full)+geom_line(aes(x=date,y=EC,color="EC"))+geom_point(aes(x=date,y=EC,color=OrigenTag))

# meteo ####
meteoobs <- read_excel("/media/usuario/32b62ac8-a81b-4630-bf83-822c71ee0cad/mdiaz/Documents/paper_facu/data/200161 BUENOS AIRES OBSERVATORIO.xlsx")
meteoobs$date=meteoobs$FECHA+meteoobs$HORA*3600
meteoobs$date_puestadefiltro=meteoobs$date-12*3600
meteoobs$dateorig=meteoobs$date
meteoobs$date=meteoobs$date_puestadefiltro
meteoobsday=timeAverage(meteoobs[,c("date","PRECIP 6HS (mm)")],avg.time = "day",statistic = "sum")
# TVmeteo=timeVariation(meteoobs,pollutant = "PRECIP 6HS (mm)")
ggplot()+geom_line(data=PMF_BA_full,aes(x=date,y=`PM2,5`))+
  geom_point(data=meteoobsday,aes(x=date,y=`PRECIP 6HS (mm)`*2, color="precip"))
TVmeteoobs=timeVariation(meteoobsday, pollutant = "PRECIP 6HS (mm)", statistic = "median")
print(TVmeteoobs$plot$month)
meteoobsmonth=timeAverage(meteoobs[,c("date","PRECIP 6HS (mm)")],avg.time = "month",statistic = "sum")
ggplot()+geom_line(data=PMF_BA_full,aes(x=date,y=`PM2,5`))+
  geom_line(data=meteoobsday,aes(x=date,y=`PRECIP 6HS (mm)`*2, color="precip"))
observaciones = merge(meteoobsday, PMF_BA_full, by="date")
observaciones = observaciones[,c("date", "PRECIP 6HS (mm)","PM2,5")]

library(ggplot2)
library(dplyr)
library(tidyr)

# Crear DataFrame combinado
observaciones <- merge(meteoobsday, PMF_BA_full, by = "date")
observaciones <- observaciones[, c("date", "PRECIP 6HS (mm)", "PM2,5")]

# Formatear y procesar el DataFrame
observaciones <- observaciones %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d")) %>% # Asegurar formato de fecha
  rename(PP = `PRECIP 6HS (mm)`) %>%                   # Renombrar columna
  # mutate(PP = replace_na(PP, 0)) %>%                   # Reemplazar NA por 0
  mutate(mes_ano = format(date, "%m-%Y")) %>%          # Crear columna mes-año
  arrange(date) %>%                                    # Ordenar por fecha
  mutate(mes_ano = factor(mes_ano, levels = unique(mes_ano))) # Orden cronológico

# Reorganizar los datos para ggplot
observaciones_long <- observaciones %>%
  pivot_longer(
    cols = c(`PP`, `PM2,5`),
    names_to = "variable",
    values_to = "value"
  )

# Crear el boxplot agrupado por mes-año
ggplot(observaciones_long, aes(x = mes_ano, y = value, fill = variable)) +
  geom_boxplot(position = position_dodge(width = 0.8),outliers = FALSE) +
  labs(
    title = "Boxplot de Precipitación y PM2.5 por Mes-Año",
    x = "Mes-Año",
    y = "Valor",
    fill = "Variable"
  ) +
  scale_fill_manual(values = c("lightblue", "lightgreen")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotar etiquetas del eje x

# observaciones = rename(observaciones, `PP`=`PRECIP 6HS (mm)`)library(ggplot2)
library(dplyr)
library(tidyr)

# Crear DataFrame combinado
observaciones <- merge(meteoobsday, PMF_BA_full, by = "date")
observaciones <- observaciones[, c("date", "PRECIP 6HS (mm)", "PM2,5")]

# Procesar y formatear datos
observaciones <- observaciones %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d")) %>% # Asegurar formato de fecha
  rename(PP = `PRECIP 6HS (mm)`) %>%                   # Renombrar columna
  mutate(PP = replace_na(PP, 0)) %>%                   # Reemplazar NA por 0
  mutate(mes_ano = format(date, "%m-%Y")) %>%          # Crear columna mes-año
  arrange(date) %>%                                    # Ordenar por fecha
  mutate(mes_ano = factor(mes_ano, levels = unique(mes_ano))) # Orden cronológico

# Reorganizar los datos para gráfico de líneas
observaciones_long <- observaciones %>%
  pivot_longer(
    cols = c(`PP`, `PM2,5`),
    names_to = "variable",
    values_to = "value"
  )

# Crear gráfico de líneas
ggplot(observaciones_long, aes(x = date, y = value, color = variable, group = variable)) +
  geom_line(size = 1) + # Líneas
  geom_point(size = 2) + # Puntos para destacar valores
  labs(
    title = "Gráfico de Líneas de Precipitación y PM2.5",
    x = "Fecha",
    y = "Valor",
    color = "Variable"
  ) +
  scale_color_manual(values = c("lightblue", "darkgreen")) + # Colores personalizados
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotar etiquetas del eje x
    legend.position = "top" # Leyenda en la parte superior
  )

library(ggplot2)
library(dplyr)
library(tidyr)

# Asegurarse de que la columna 'date' esté en formato Date y extraer el mes
observaciones <- observaciones %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d")) %>%
  mutate(month = factor(format(date, "%m"), levels = sprintf("%02d", 1:12)))

# Renombrar la columna y reemplazar NA por 0
observaciones <- observaciones %>%
  rename(PP = `PRECIP 6HS (mm)`) %>%
  mutate(PP = replace_na(PP, 0))

# Reorganizar los datos para ggplot
observaciones_long <- observaciones %>%
  pivot_longer(
    cols = c(`PP`, `PM2,5`),
    names_to = "variable",
    values_to = "value"
  )

# Crear el boxplot con las dos variables
ggplot(observaciones_long, aes(x = month, y = value, fill = variable)) +
  geom_boxplot(position = position_dodge(width = 0.8),outliers = FALSE) +
  labs(
    title = "Boxplot de Precipitación y PM2.5 por Mes",
    x = "Mes",
    y = "Valor",
    fill = "Variable"
  ) +
  scale_fill_manual(values = c("lightblue", "lightgreen")) +
  theme_minimal()


library(ggplot2)
library(lubridate)
library(dplyr)

# Crear una columna con el mes
# Asegurarse de que 'date' esté en formato Date
observaciones <- observaciones %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d")) %>%
  mutate(month = factor(format(date, "%m"), levels = sprintf("%02d", 1:12)))

# Crear el boxplot
ggplot(observaciones, aes(x = month, y = `PM2,5`)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  labs(
    title = "Boxplot de PM2.5 por Mes",
    x = "Mes (número)",
    y = "PM2.5 (µg/m³)"
  ) +
  theme_minimal()
ggplot()+
  geom_line(data=meteoobsmonth,aes(x=date,y=`PRECIP 6HS (mm)`, color="precip"))



monthPM= timeAverage(PMF_BA_full[,c("date","PM2,5")], avg.time = "month", statistic="median")
ggplot()+geom_line(data=monthPM,aes(x=date,y=`PM2,5`, color="PMmensual"))+
  geom_line(data=meteoobsmonth,aes(x=date,y=`PRECIP 6HS (mm)`/10, color="acumuladamensual"))

observaciones <- observaciones %>%
  mutate(precip_rango = case_when(
    `PRECIP 6HS (mm)` == 0 ~ "Sin precipitación",
    `PRECIP 6HS (mm)` > 0 & `PRECIP 6HS (mm)` <= 10 ~ "Baja",
    `PRECIP 6HS (mm)` > 10 & `PRECIP 6HS (mm)` <= 50 ~ "Moderada",
    `PRECIP 6HS (mm)` > 50 ~ "Alta"
  ))

ggplot(observaciones, aes(x = month, y = `PM2,5`, fill = precip_rango)) +
  geom_boxplot() +
  facet_wrap(~ precip_rango, ncol = 2) +
  labs(
    title = "PM2.5 por Rango de Precipitación y Mes",
    x = "Mes",
    y = "PM2.5 (µg/m³)",
    fill = "Rango de Precipitación"
  ) +
  theme_minimal()



library(ggplot2)
library(dplyr)

# Calcular precipitación promedio por mes
precip_avg <- observaciones %>%
  group_by(month) %>%
  summarise(precip_promedio = mean(`PRECIP 6HS (mm)`, na.rm = TRUE))

# Combinar datos
observaciones <- observaciones %>%
  left_join(precip_avg, by = "month")

# Crear el boxplot con color por precipitación promedio
ggplot(observaciones, aes(x = month, y = `PM2,5`, fill = precip_promedio)) +
  geom_boxplot(color = "black") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title = "Relación entre PM2.5 y Precipitación por Mes",
    x = "Mes",
    y = "PM2.5 (µg/m³)",
    fill = "Precipitación Promedio (mm)"
  ) +
  theme_minimal()

library(ggplot2)
library(patchwork)

# Gráfico de precipitación promedio
g1 <- ggplot(precip_avg, aes(x = month, y = precip_promedio)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(
    title = "Precipitación Promedio por Mes",
    x = "Mes",
    y = "Precipitación (mm)"
  ) +
  theme_minimal()

# Gráfico de boxplot para PM2.5
g2 <- ggplot(observaciones, aes(x = month, y = `PM2,5`)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  labs(
    title = "PM2.5 por Mes",
    x = "Mes",
    y = "PM2.5 (µg/m³)"
  ) +
  theme_minimal()

# Combinar los gráficos
g1 / g2
# otros relaciones ####

ggplot(PMF_BA_full[PMF_BA_full$Ni!=0,])+geom_point(aes(x=Ni,y=V))+
  geom_abline(aes(slope=2.4, intercept=0), linetype="dashed")+
  geom_abline(aes(slope=3, intercept=0, color="oil combustion from ship engines"), linetype="dashed")+
  geom_abline(aes(slope=1.49, intercept=0, color="road traffic"), linetype="dashed")+
  geom_abline(aes(slope=0.01, intercept=0, color="prueba"), linetype="dashed")

ggplot(PMF_BA_full[PMF_BA_full$Cd>0.00005 & PMF_BA_full$Sb>0,])+geom_point(aes(x=Cd,y=Sb))+
  geom_abline(aes(slope=5, intercept=0, color="road brake pad wear "), linetype="dashed")+
  geom_abline(aes(slope=8, intercept=0, color="prueba"), linetype="dashed")

ggplot(PMF_BA_full[PMF_BA_full$Cd<0.005 & PMF_BA_full$Cd>0.00001,],aes(x=Cd,y=Sb))+geom_point()+
  stat_summary(fun.data= mean_cl_normal) + 
  geom_smooth(method='lm',se = F)+
  coefs <- coef(lm(Sb~Cd, data = PMF_BA_full_so[PMF_BA_full_so$Cd>0.00001,]))


# PMF
X2C5f_base <- read_excel("~/Documents/paper_facu/data/2C5f_base.xlsx", 
                         sheet = "Contributions", range = "B108:H208")
X2C5f_base <- dplyr::rename(X2C5f_base, date = '...1')
TV=timeVariation(X2C5f_base,pollutant = "Factor 2")


library(readr)
blhERA5 <- read_csv("/media/usuario/32b62ac8-a81b-4630-bf83-822c71ee0cad/mdiaz/Documents/paper_facu/data/DatosCNEAgases/blh20192020CNEACent.txt")
blhERA5$date=blhERA5$date-3600*12
blhERA5day <- timeAverage(blhERA5,avg.time = "day")
mergeblhPm2.5=merge(blhERA5day,PMF_BA_full, by="date")
mergeblhPm2.5=merge(mergeblhPm2.5,BA_events_testM, by="date")
ggplot(mergeblhPm2.5, aes(x=`PM2,5`,y=blh))+geom_point()
mergeblhPm2.5$wspd=sqrt(mergeblhPm2.5$u10^2+mergeblhPm2.5$v10^2)
mergeblhPm2.5$ventcoeff=mergeblhPm2.5$blh*mergeblhPm2.5$wspd
ggplot(mergeblhPm2.5, aes(x=`PM2,5`,y=ventcoeff,color=Event_F))+geom_point()
ggplot(mergeblhPm2.5, aes(x=`PM2,5`,y=blh,color=Event_F))+geom_point()
ggplot(mergeblhPm2.5, aes(x=`PM2,5`,y=wspd,color=Event_F))+geom_point()
