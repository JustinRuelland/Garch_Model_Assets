#--------------Imports des packages--------------------
rm(list=ls())

install.packages("ggplot2")
install.packages("dplyr")

library (tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate)

#---------------Imports des autres fonctions---------------------
source(file = "./data_preparation.R",local= TRUE)
source(file= "./condition_stationnarite.R",local=TRUE)


#---------------Import et transformation des données--------------
data = read.csv("./^GDAXI.csv") #mettre le fichier de votre choix (provenant de Yahoo finance)

data <- transform_csv(data)

#Premiers graphiques - statistiques descriptives
plot_series_temp(data)
autocorrelations(data)


#--------------
cond_statio(rnorm,0,1)







#------------------appel d'une fonction d'un autre fichier.R---------------------------

#mettez l'url de votre pc vers le fichier (mettre des / et pas des \ dans l'url)
source(file= "./condition_stationnarite.R",local=TRUE)
#on peut executer une fonction de ce fichier importé, par exemple "cond_statio":
cond_statio(rnorm,0,1)
