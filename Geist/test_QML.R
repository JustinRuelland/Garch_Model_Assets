# Test du QML pour application dans puissance_test_chgtGARCH.R - 

library("dplyr")
library("purrr")
library("plyr")
library("ggplot2")

source(file="./simulation_series.R",local=TRUE)
source(file= "./QML_Variance.R",local=TRUE)
source(file= "./prevision.R",local=TRUE)

# Options de la fonction
cut_chgt=0.01
n_path=1
n=100000

# Fonctions chgt GARCH
niveau_test = 0.05
theta1 = c(0.0001,0.12,0.85)

# Génération des données
alpha_beta = c(0,0)

rendements = simulation_rendements_avec_changement_GARCH(n,theta1,unlist(c(0.0001,c(alpha_beta))),cut_chgt)

# Appel de QML (qui dans la fonction de base se fait dans func_backtest)
theta = QML(rendements**2)

print(theta)
