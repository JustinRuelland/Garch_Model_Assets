install.packages("signal")
install.packages("tidyverse") #contient notamment ggplot2, dplyr
install.packages("lubridate")

library(signal)
library (tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate)


library(zoo)
library(roll)
library(ggplot2)
library(dplyr)
library(pracma)
#install.packages("dm.test")
#install.packages('forecast', dependencies = TRUE)
library('dm.test')
#install.packages("ggplot2")          # Install ggplot2 package
library("ggplot2")   

#---------------------- Mise en place des données ------------------------
source(file = "C:/Users/thoma/Documents/GitHub/Garch_Model_Assets/data_preparation.R",local= TRUE)
data = read.csv("C:/Users/thoma/Documents/GitHub/Garch_Model_Assets/CAC40_15_19.csv") #fichier csv de Yahoo finance

transform_csv <- function(data){
  data = select(data, "Open","Date")
  data = rename(data, c("Prix"="Open"))
  data$Prix = as.numeric(data$Prix)
  data$Date = as.Date(data$Date)
  
  data = mutate(data, rendement = log(data$Prix/lag(data$Prix)))
  data = mutate(data, rendement2 = rendement**2)
  return(data[-1,])
}

data <- transform_csv(data)



#---------------------- Calcul du sigma par QML ------------------------

#on part du principe qu'on a observe une série financière
#dont les rendements au carré sont eps2

QML <-function(eps2){
  n = length(eps2) #longueur de la série financière
  
  #fonction de vraisemblance à optimiser
  f_opt <- function(theta_opt){
    sigmas2_QML = simu_sigma2(eps2,theta_opt) #simulation des sigmas à partir des eps2
    #on retire les 20 premières données (négligeables cf notes)
    return(sum(log(sigmas2_QML[25:n])+(eps2[25:n]/sigmas2_QML[25:n]))) } #log vraisemblance
  
  theta_init = c(0.0001,0.12,0.8) #valeur initiale
  
  #contraintes
  ui <- cbind(c(1,-1,0,0,0,0),c(0,0,1,-1,0,0),c(0,0,0,0,1,-1))
  ci <- c(10**(-9), -1, 0, -3, 0, -0.99)
  return(constrOptim(theta=theta_init,f = f_opt,ci=ci,ui=ui,gr=NULL)$par)} #opti sous contraintes linéaires



#simulation des sigmas**2 à partir des epsilon**2 connus
simu_sigma2 <- function(eps2,theta){
  n = length(eps2)
  
  sigmas2 = double(n) #création d'un vecteur de taille qui ne va pas varier par la suite
  sigmas2[1] = c(mean(eps2))
  
  for(i in 2:n){sigmas2[i]=theta[1]+theta[2]*eps2[i-1]+theta[3]*sigmas2[i-1]}
  
  return(sigmas2)
}


#calcul récurssif pour un sigma2(t) précis
func_sigma2 <- function(t,sigma2_init,eps2,theta){
  res <- sigma2_init
  if(t>1){
    res<-theta[1] + theta[2]*eps2[t-1] + theta[3]*func_sigma2(t-1,sigma2_init,eps2,theta)}
  return(res)}

#---------------------- Prediction des rendements par GARCH à horizon 1  ------------------------


pred_horizon_1 <- function(cut, eps2){
  n = length(eps2)
  df <- data.frame(Col1 = double())
  theta = QML(eps2[1:cut])
  sigma2 = simu_sigma2(eps2[1:cut], theta)
  t = cut+1
  for(i in t:n) {
    sigma2_t = func_sigma2(i, sigma2[1], eps2[1:i-1], theta)
    df[i,1] <-sigma2_t 
  }
  return(data.frame(df))
}

pred_horizon_1_cut <- function(cut, eps2){
  n = length(eps2)
  df <- data.frame(Col1 = double())
  t = cut+1
  
  for(i in t:n) {
    borne = i-1
    theta = QML(eps2[1:borne])
    sigma2 = simu_sigma2(eps2[1:borne], theta)
    sigma2_t = func_sigma2(i, sigma2[1], eps2[1:borne], theta)
    df[i,1] <-sigma2_t 
  }
  return(data.frame(df))
}


cut = 500
eps2 = data[,4]

pred = data.frame(pred_horizon_1(cut, eps2))
pred = data.frame(pred[501:1021,1])

pred_cut = data.frame(pred_horizon_1_cut(cut, eps2))
pred_cut = data.frame(pred_cut[501:1021,1])

eps2 = data.frame(data[501:1021,4])


#---------------------- Prédiction des rendements par méthode naïve à horizon 1 ------------------------
bb_exp<-function(data, window, m){
  
  sma_na = data.frame(x= c())
  upper_na = data.frame(x= c())
  lower_na = data.frame(x= c())
  for(i in 1:window){
    sma_na[i,1] = NA
    upper_na[i,1] = NA
    lower_na[i,1] = NA
  }
  
  
  sma = movavg(data, n=window, type='e')[window:length(movavg(data, n=window, type='e'))]
  std = roll_sd(data, window )[window:length(roll_sd(data, window ))]
  upper_bb = sma + std * m
  lower_bb = sma - std * m
  
  
  for(i in 1:length(sma)){
    sma_na[window+i,1] = sma[i]
    upper_na[window+i,1] = upper_bb[i]
    lower_na[window+i,1] = lower_bb[i]
  }
  
  
  vecteur1 = data.frame(data)
  vecteur1[length(data)+1,1] = NA
  
  tableau = data.frame(x = vecteur1, y = upper_na, z = lower_na, a = sma_na) 
  return(tableau)
}

bb_data = bb_exp(data[,4],4, 2)
bb_data = data.frame(bb_data[501:1021,4])



#---------------------- Comparaison des erreurs ------------------------
data2 <- data.frame(x = 1:length(eps2[,1]),
                   y1 = eps2[,1],
                   y2 = pred[,1],
                   y3 = pred_cut[,1],
                   y4 = bb_data[,1])


ggp1 <- ggplot(data2, aes(x)) +       # Create ggplot2 plot
  geom_line(aes(y = y1), color = "red") +
  geom_line(aes(y = y2), color = "blue") +
  geom_line(aes(y = y3), color = "green") +
  geom_line(aes(y = y4), color = "black")
ggp1                  

#---------------------- Test DM  ------------------------
#---------------------- For alternative="greater", the alternative hypothesis is that method 2 is more accurate than method 1. ------------------------


forecast::dm.test(pred[,1]-eps2[,1], pred_cut[,1]-eps2[,1], alternative = c("two.sided"), h = 1, power = 2)
forecast::dm.test(pred[,1]-eps2[,1], pred_cut[,1]-eps2[,1], alternative = c("greater"), h = 1, power = 2)

forecast::dm.test(bb_data[,1]-eps2[,1], pred[,1]-eps2[,1], alternative = c("two.sided"), h = 1, power = 2)
forecast::dm.test(bb_data[,1]-eps2[,1], pred_cut[,1]-eps2[,1], alternative = c("two.sided"), h = 1, power = 2)
forecast::dm.test(bb_data[,1]-eps2[,1], pred_cut[,1]-eps2[,1], alternative = c("greater"), h = 1, power = 2)

