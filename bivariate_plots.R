library(openair)
library(ggplot2)
library(lubridate)
setwd("/home/pablo/doctorado/paper_facu")

data <- read.csv("~/doctorado/paper_facu/data_every_hour_obs.csv")
data$date <- as.POSIXct(data$date, tz='UTC')
print(typeof(data$date))

keys = names(data)[2:45]

data_period2 = data[which(data$date >= as.POSIXct("2019-05-22", format="%Y-%m-%d") & data$date <= as.POSIXct("2019-06-04", format="%Y-%m-%d")), ]
print(data_period2)
windRose(data_period2)#, poll = 'Na.sol', stati = 'cpf', percentile=c(0, 100),  main=paste('Na sol','cpf'))

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


polarPlot(data)


data_event = data[which(data$Event%in%c('S', 'SP', 'SN')),]
max_perc = 25
min_perc = 0
for (title in keys) {
  print(title)
  png(paste('images/CPF/CPF_Q1_',title,'.png', sep=''))
  polarPlot(data, poll = title, stati = 'cpf', percentile = c(min_perc,max_perc), main=paste(title,'Q1'))
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

date <- as.POSIXct("2019-17-7", format="%Y-%m-%d")
data_period1 = data[which(data$date < as.POSIXct("2019-07-17", format="%Y-%m-%d")),]
data_period2 = data[which(data$date >= as.POSIXct("2019-07-17", format="%Y-%m-%d") & data$date < as.POSIXct("2019-11-14", format="%Y-%m-%d")), ]
print(data_period2)

polarPlot()

print("Periods")
for (title in keys) {
  print(paste(title, "Period 1"))
  png(paste("images/CPF/Na_Period1/CPF_Period1_", title, "_Q4.png", sep=""))
  polarPlot(data_period1, poll = title, stati = 'cpf', percentile = c(75, 100), main=paste(title, 'before Jul17 Q4'))
  dev.off()
  print('done Q4')
  if (title=='Co'){
    png(paste("images/CPF/Na_Period1/CPF_Period1_", title, "_Q1.png", sep=""))
    polarPlot(data_period1, poll = title, stati = 'cpf', percentile = c(0, 25), main=paste(title, 'before Jul17 Q1'), limits=c(0, 8.441e-05))
    dev.off()
  } else {
    png(paste("images/CPF/Na_Period1/CPF_Period1_", title, "_Q1.png", sep=""))
    polarPlot(data_period1, poll = title, stati = 'cpf', percentile = c(0, 25), main=paste(title, 'before Jul17 Q1'))
    dev.off()
  }

  print('done Q1')
  print(paste(title, "Period 2"))
  png(paste("images/CPF/Na_Period2/CPF_Period2", title, "_Q4.png", sep=""))
  polarPlot(data_period2, poll = title, stati = 'cpf', percentile = c(75, 100), main=paste(title, 'Jul17 to Nov14 Q4'))
  dev.off()
  print('done Q4')
  if (title=='Ca'){
    png(paste("images/CPF/Na_Period2/CPF_Period2", title, "_Q1.png", sep=""))
    polarPlot(data_period2, poll = title, stati = 'cpf', percentile = c(0, 25), main=paste(title, 'Jul17 to Nov14 Q1'), limits=c(0, 0.025))
    dev.off()
  } else if (title=='Co'){
    png(paste("images/CPF/Na_Period2/CPF_Period2", title, "_Q1.png", sep=""))
    polarPlot(data_period2, poll = title, stati = 'cpf', percentile = c(0, 25), main=paste(title, 'Jul17 to Nov14 Q1'), limits=c(0, 8.441e-05))
    dev.off()
  } else {
    png(paste("images/CPF/Na_Period2/CPF_Period2", title, "_Q1.png", sep=""))
    polarPlot(data_period2, poll = title, stati = 'cpf', percentile = c(0, 25), main=paste(title, 'Jul17 to Nov14 Q1'))
    dev.off()
  }
    print('done Q1')
}

print(data_period2$Ca[which(data_period2$Ca<=0.0250)])
