#--------------Imports des packages--------------------
rm(list=ls())

install.packages("ggplot2")
install.packages("dplyr")

library (tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate)


#---------------Import et transformation des données--------------
source(file = "./data_preparation.R",local= TRUE)

data = read.csv("./^GDAXI.csv") #mettre le fichier de votre choix (provenant de Yahoo finance)

data <- transform_csv(data)

#Premiers graphiques - statistiques descriptives
plot_series_temp(data)
autocorrelations(data)


#--------------- Conditions de stationnarité ---------------------
source(file= "./condition_stationnarite.R",local=TRUE)
condition_stationnarite(rnorm)

runif_normalisee <-function(n){return(runif(n,-sqrt(3),sqrt(3)))}
condition_stationnarite(runif_normalisee)

rt_8_normalisee <-function(n){return(rt(n,8)/sqrt(8/(8-2)))} # Loi de Student à 8 degrés de liberté
condition_stationnarite(rt_8_normalisee)

#A faire #Loi de Mises



#------------------appel d'une fonction d'un autre fichier.R---------------------------

#mettez l'url de votre pc vers le fichier (mettre des / et pas des \ dans l'url)
source(file= "./condition_stationnarite.R",local=TRUE)
#on peut executer une fonction de ce fichier importé, par exemple "cond_statio":
condition_stationnarite(rnorm)