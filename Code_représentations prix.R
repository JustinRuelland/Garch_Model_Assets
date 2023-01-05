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

## Prix & rendements
prix = ggplot(data = DAX) + geom_line(aes(x = Date,y = Open))
plot_rdt = ggplot(data = DAX) + geom_line(aes(x = Date,y = rendement))

plot(prix)
plot(plot_rdt)

## Auto-corrélations
acf(na.omit(DAX$rendement))
acf(na.omit((DAX$rend_car)))

### Simulations

## Représentations des (alpha,beta) vérifiant la condition de stationnarité
n = 10**5


Condition_stricte_statio <-function(x){
res <- c()
x<-x**2

Abscisse_alpha=c(0:199)/100
for(alpha in Abscisse_alpha){
  
  for(b in c(1:200)){
    beta = (b-1)/100
    
    esp = mean(log(alpha*x + beta))
    
    if (esp>=0){
      res <- c(res,beta)
      break
    }
    
  }
}
Beta_values = res

# Création du dataframe pour le mettre dans le ggplot

Alpha_Beta = cbind(Abscisse_alpha,Beta_values)
res_avec_abscisse = as.data.frame(Alpha_Beta)

# ggplot
Graph_condi_satio = ggplot(data = res_avec_abscisse) + aes(ymax = 1,ymin=0)+ geom_line(aes(x = Abscisse_alpha,y = Beta_values))
print(Graph_condi_satio)

return(1)
}


# Loi normale
x = rnorm(n)
Condition_stricte_statio(x)

# Loi de Cauchy
x = rcauchy(n)
Condition_stricte_statio(x)

# Loi uniformes
x = runif(n,min=-1,1)*sqrt(3)
Condition_stricte_statio(x)

# Loi Mises
x = (rexp(n,2)-1/2)*2
Condition_stricte_statio(x)
