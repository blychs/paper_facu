library(openair)
library(ggplot2)
library(lubridate)
library(readxl)
BA_events_testM <- read_excel("BA_events_testM.xlsx",
                              col_types = c("date", "numeric", "numeric",
                                            "skip", "text"))
BA_events_testM$date <- as.POSIXct(BA_events_testM$date, tz='UTC')
BA_events_testM$Event_M <- as.factor(BA_events_testM$Event_M)

factor(BA_events_testM$Event_M)

PMF_BA_full <- read_excel("PMF_BA_full.xlsx",
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
PMF_BA_full$OC_EC = PMF_BA_full$OC/PMF_BA_full$EC
PMF_BA_full$OC_K=PMF_BA_full$OC/PMF_BA_full$K
PMF_BA_full_so=PMF_BA_full[PMF_BA_full$date<=as.POSIXct('2019-05-23') | PMF_BA_full$date>=as.POSIXct('2019-06-02'),]
corPlot(PMF_BA_full_so, dendrogram = TRUE,method = "pearson", main= "R pearson, sin outliers")
corPlot(PMF_BA_full, dendrogram = TRUE,method = "pearson", main ="R pearson, matriz completa")

corPlot(PMF_BA_full, dendrogram = TRUE,method = "spearman", main ="R spearman, matriz completa")
corPlot(PMF_BA_full_so, dendrogram = TRUE,method = "spearman", main ="R spearman, sin outliers")

# Define event levels
event_levels <- c("SN", "SL", "S")

# Count the number of events
event_count <- sum(table(BA_events_testM$Event_M[month(BA_events_testM$date)>=3&month(BA_events_testM$date)<=5])[event_levels])
print(event_count)
event_count <- sum(table(BA_events_testM$Event_M[month(BA_events_testM$date)>=6&month(BA_events_testM$date)<=8])[event_levels])
print(event_count)
event_count <- sum(table(BA_events_testM$Event_M[month(BA_events_testM$date)>=9&month(BA_events_testM$date)<=11])[event_levels])
print(event_count)
event_count <- sum(table(BA_events_testM$Event_M[month(BA_events_testM$date)==12|month(BA_events_testM$date)<=2])[event_levels])
print(event_count)

