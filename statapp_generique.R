#------script générique pour le projet statapp--------------------

install.packages("ggplot2")
library(ggplot2)
install.packages("dplyr")
library(dplyr)

#---------------importation et transformation des données--------------

data = read.csv("Ensae/2A/Statapp/donnees/data_cac40.csv") #mettre le fichier de votre choix (provenant de Yahoo finance)
data = select(data, "Open")
data = rename(data, c("Prix"="Open"))
data = mutate(data, eps = log(data$Prix/lag(data$Prix)))#rendements (eps)
data = mutate(data, eps2 = eps**2) #rendements au carré (eps2)
temp = seq_len(dim(data)[1]) #temps


#-----------------premiers graphiques ----------------
plot(data$Prix, type = "l", main = "Prix de l'actif")
plot(data$eps, type="l", col="blue", main = "Rendements de l'actif")
plot(data$eps2, type="l", col="black", main = "Rendements au carré de l'actif")

acf(data$eps,type='correlation', na.action=na.pass, plot=TRUE)
acf(data$eps2,type='correlation', na.action=na.pass, plot=TRUE)

#------------------appel à une fonction d'un autre fichier.R---------------------------

source(file= "C:/Users/maeld/OneDrive/Documents/GitHub/Garch_Model_Assets/fonctions_statio.R",local=TRUE)
cond_statio(rnorm,0,1)
