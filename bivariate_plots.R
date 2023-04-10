library(openair)
library(ggplot2)
library(lubridate)
setwd("/home/pablo/doctorado/paper_facu")

data <- read.csv("~/doctorado/paper_facu/data_every_hour_obs.csv")
data$date <- as.POSIXct(data$date, tz='UTC')
print(typeof(data$date))

keys = names(data)[2:45]

#data[date(data$date) > as.POSIXct('2019-05-23') & date(data$date) < as.POSIXct('2019-06-04') ,-1] = NA
#data[date(data$date) == as.POSIXct(c('2020-02-21', '2020-02-22')),] = NA
max_perc = 100
min_perc = 75
for (title in keys) {
print(title)
png(paste('images/median/bv_median_',title,'.png', sep=''))
polarPlot(data, poll = title, stati = 'median',  main=paste(title,'median'))
dev.off()
}


data_event = data[which(data$Event%in%c('S', 'SP', 'SN')),]
max_perc = 100
min_perc = 75
for (title in keys) {
  print(title)
  png(paste('images/CPF/CPF_Event_P75_',title,'.png', sep=''))
  polarPlot(data_event, poll = title, stati = 'cpf', percentile = c(min_perc,max_perc), main=paste(title,'P75'))
  dev.off()
}

data_no_event = data[which(!data$Event%in%c('S', 'SP', 'SN')),]
max_perc = 100
min_perc = 75
for (title in keys) {
  print(title)
  png(paste('images/CPF/CPF_No_Event_P75_',title,'.png', sep=''))
  polarPlot(data_no_event, poll = title, stati = 'cpf', percentile = c(min_perc,max_perc), main=paste(title,'P75'))
  dev.off()
}
