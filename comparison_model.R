rm(list=ls())

install.packages("forecast")
library(forecast)

install.packages("signal")
install.packages("tidyverse") #contient notamment ggplot2, dplyr
install.packages("lubridate")

library(signal)
library (tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate)

source(file= "./QML_Variance.R",local=TRUE)
source(file = "./data_preparation.R",local= TRUE)
data = read.csv("./CAC40_15_19.csv") #fichier csv de Yahoo finance
data <- transform_csv(data)

test_mariano <- function(pred1,pred2,val,loss,hor){
  e1 = pred1-val
  e2 = pred2-val
  test =  dm.test(pred1-val, pred2-val, 
              alternative = c("two.sided"), h = hor, power = 2)
  
  res=data.frame(X=c(test.p.value,test.statistic))
  Y=c('p_value','stat_test')
  res$Y=Y[res$X]
  
  return(res)}

pred_h1_garch <- function(eps2,cut){
  n = length(eps2)
  theta =  QML(eps2[1:floor(n*cut)])
  pred = c()
  pred[1] = theta[1]/(1-theta[2]-theta[3])
  for(i in 2:n-floor(n*cut)-1){
    print(eps2[floor(n*cut)+i-1])
    pred[i] = theta[1]+theta[2]*eps2[floor(n*cut)+i-1]+theta[3]*pred[i-1]}
  return(pred)}

pred_h1_garch(data$rendement2,0.8)
