rm(list=ls())

install.packages("signal")
install.packages("tidyverse") #contient notamment ggplot2, dplyr
install.packages("lubridate")

library(signal)
library (tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate)


library(zoo)
install.packages("forecast")
library(forecast)


source(file= "./QML_Variance.R",local=TRUE)
source(file = "./data_preparation.R",local= TRUE)
data = read.csv("./CAC40_15_19.csv") #fichier csv de Yahoo finance
data <- transform_csv(data)

test_mariano <- function(pred1,pred2,val,hor){
  e1 = pred1-val
  e2 = pred2-val
  test =  dm.test(pred1-val, pred2-val, 
              alternative = c("two.sided"), h = hor, power = 2)
  
  return(test)}

pred_h1_garch <- function(eps2,cut){
  
  n = length(eps2)
  n_cut = floor(n*cut)
  
  theta =  QML(eps2[1:n_cut])
  pred = double(n-n_cut)
  init = theta[1]/(1-theta[2]-theta[3])
  pred[1] = func_sigma2(n_cut,init,eps2[1:n_cut],theta)
  
  l = n-n_cut
  for(i in 2:l){
    pred[i] = theta[1]+theta[2]*eps2[n_cut+i-2]+theta[3]*pred[i-1]}
  return(pred)}


rolling_av <- function(eps2,cut){
  
  n = length(eps2)
  n_cut = floor(n*cut)
  res_temp = rollmean(x=c(rep(NA,9),eps2),k = 10, align='right',na.rm = TRUE)
  print(res_temp)
  s = n_cut+1
  return(res_temp[s:length(res_temp)])}



l = floor(0.8*length(data$rendement2))
yy_garch = pred_h1_garch(data$rendement2,0.8)
plot(data$rendement2[l:length(data$rendement2)],type='l')
lines(c(1:205),yy_garch, col='red')

plot(data$rendement2[l:length(data$rendement2)],type='l')
yy_roll = rolling_av(data$rendement2,0.8)
lines(c(1:205),yy_roll, col='red')

test_mariano(yy_garch,yy_roll,data$rendement2[l:length(data$rendement2)],hor=1)
