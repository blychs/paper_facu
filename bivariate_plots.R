library(openair)
library(ggplot2)
setwd("/home/pablo/doctorado/paper_facu")

data <- read.csv("~/doctorado/paper_facu/data_every_hour_duplicate.csv")
print(names(data))

keys = names(data)[2:45]

max_perc = 100
min_perc = 75
for (title in keys) {
print(title)
png(paste('images/CPF/CPF_P75_',title,'.png', sep=''))
polarPlot(data, poll = title, stati = 'cpf', percentile = c(min_perc,max_perc), main=paste(title,'P75'))
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
