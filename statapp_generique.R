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

#------------------appel d'une fonction d'un autre fichier.R---------------------------

#mettez l'url de votre pc vers le fichier (mettre des / et pas des \ dans l'url)
source(file= "./fonctions_statio.R",local=TRUE)
#on peut executer une fonction de ce fichier importé, par exemple "cond_statio":
cond_statio(rnorm,0,1)
