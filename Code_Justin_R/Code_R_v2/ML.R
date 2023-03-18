rm(list=ls())
setwd(dir = "/Users/justinr/Documents/2A/Stat App/Github/Garch_Model_Assets/Code_Justin_R/Code_R_v2")
library(ggplot2)
library(dplyr)

source(file="./../../simulation_series.R",local=TRUE)
source(file="./../../QML_Variance.R",local=TRUE)
source(file="./../../data_preparation.R",local=TRUE)

# Regression linéaire 

# Import series financieres reelles
data_real <- read.csv("./../../^GDAXI.csv") #fichier csv de Yahoo finance
data_real <- transform_csv(data_real)

eps_square <- data_real$rendement2

# Dropping na
eps_square <- eps_square[!is.na(eps_square)]

# Visualisation des correlograms pour determiner les variables d'interets
autocorrelations_korrigieren <-function(data){
  #acf(data$rendement,type='correlation', na.action=na.pass, lag.max = 0,plot=TRUE)
  acf(data$rendement2,type='correlation', na.action=na.pass, lag.max = 100,plot=TRUE)
}


autocorrelations_korrigieren(data_real)
# On identifie une baisse de correlation au dela de 40

# Creation du database

colname_erstellung <- function(lag.max){
  naming <- as.character(1:lag.max)
  for(i in 1:lag.max){
    naming[i] <- paste(naming[i], "lag")
  }
  return(naming)
}

database_erstellung <- function(serie, lag.max){
  taille <- length(serie)
  nb_ligne <- taille - lag.max
  data_matrix_X <- matrix(NA, ncol = lag.max, nrow = nb_ligne)
  data_matrix_Y <- matrix(NA, ncol = 1, nrow = nb_ligne)
  debut <- lag.max + 1
  for(i in 1:nb_ligne){
    current <- i + lag.max
    data_matrix_Y[i] <- serie[current]
    first <- current - 1
    last <- current - lag.max
    data_matrix_X[i,] <- serie[first:last]
  }
  colnames(data_matrix_Y) <- "Y"
  colnames(data_matrix_X) <- colname_erstellung(lag.max)
  
  df <- as.data.frame(cbind(data_matrix_Y, data_matrix_X))
  return(df)
}

database_erstellung_prediction <- function(serie, lag.max){
  taille <- length(serie)
}

# Enregistrement design matrix
df_design <- database_erstellung(eps_square, 40)

# OLS

Y <- as.matrix(df_design$Y)
X <- as.matrix(df_design[colnames(df_design) != "Y"])
OLS_spurious <- lm(Y ~ X)

summary(OLS_spurious) 


# OLS prediction 

OLS_prediction <- function(train_set, evaluate_set, lag.max){
  df_design <- database_erstellung(train_set, lag.max)
  Y <- as.matrix(df_design$Y)
  X <- as.matrix(df_design[colnames(df_design) != "Y"])
  OLS_spurious <- lm(Y ~ X)
  
}


