rm(list=ls())

library(signal)
library (tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate)

#install.packages("forecast")
#install.packages("zoo")
#install.packages("roll")

library(zoo)
library(roll)
library(pracma)
library(forecast)
#---------------------- Mise en place des données -----------------------
data = read.csv("../CAC40_15_19.csv")
source(file="../data_preparation.R")
source(file="../QML_Variance.R")

data <- transform_csv(data)


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

