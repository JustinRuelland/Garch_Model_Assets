rm(list=ls())

install.packages("tidyverse")
setwd("C:\\Users\\louis\\Garch_Model_Assets")
DAX = read.csv("^GDAXI.csv",sep=",")

library (tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate)

### Créations des variables manipulées

DAX$Date = as.Date(DAX$Date)
DAX$Open = as.numeric(DAX$Open)
DAX = select(DAX,-"Close",-"High",-"Low",-"Close",-"Adj.Close",-"Volume")
DAX = mutate(DAX, rendement = log(DAX$Open/lag(DAX$Open))) #rendement
DAX = mutate(DAX, rend_car = rendement**2) #rendement au carré


prix = ggplot(data = DAX) + geom_line(aes(x = Date,y = Open))
plot_rdt = ggplot(data = DAX) + geom_line(aes(x = Date,y = rendement))

plot(prix)
plot(plot_rdt)

acf(na.omit(DAX$rendement))
acf(na.omit((DAX$rend_car)))

### Simulation
n = 10**5
x = rnorm(n)
x <- x**2
res <- c()

for(a in c(1:200)){
  alpha = (a-1)/100
  
  for(b in c(1:200)){
    
    beta = (b-1)/100
    mu = mean(log(alpha*x + beta))
    
    if (mu>=0){
      res <- c(res,mu)
      break
    }
  }
}


plot(res,type="l")
