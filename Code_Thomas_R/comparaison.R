
#---------------------- Installation des packages et librairies ------------------------

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
install.packages("dm.test")
install.packages('forecast', dependencies = TRUE)
library('dm.test')

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


simu_sigma2 <- function(eps2,theta){
  n = length(eps2)
  
  sigmas2 = double(n) #création d'un vecteur de taille qui ne va pas varier par la suite
  sigmas2[1] = c(mean(eps2))
  
  for(i in 2:n){sigmas2[i]=theta[1]+theta[2]*eps2[i-1]+theta[3]*sigmas2[i-1]}
  
  return(sigmas2)
}




#---------------------- Prediction des rendements par GARCH à horizon 1  ------------------------
pred_horizon_1 <- function(cut, data, loi_eta =rnorm){
  n = length(data[,4])- cut
  df <- data.frame(Col1 = double())
  
  for(i in 1:n) {
    theta_hat = QML(data[,4][1:cut+i-1])
    
    sigma_init = sqrt(theta_hat[1]/(1-theta_hat[2]-theta_hat[3]))
    sigmas2 = c(sigma_init**2)
    etas = rnorm(1)
    etas2 = etas**2
    
    epsilons2 = etas2*sigmas2
    epsilons = sign(etas)*sqrt(epsilons2)
    df[i+cut,1] <- epsilons
  }
  return(data.frame(df))
}

cut = 800
n = length(data[,4]) - cut
pred = pred_horizon_1(800, data, loi_eta =rnorm)
a = 801
pred = data.frame(pred[801:1021,1])
epsi = data.frame(data[801:1021,3])



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

bb_data = bb_exp(data[,3],100, 2)
bb_data = data.frame(bb_data[801:1021,4])



#---------------------- Comparaison des erreurs ------------------------

x  <- c(1:221)

#plot the first data series using plot()
plot(x, pred[,1]-epsi[,1], type="o", col="blue", pch=".", ylab="y", lty=1)


#add second data series to the same chart using points() and lines()
#points(x, epsi[,1], col="red", pch=".")
#lines(x, epsi[,1], col="red",lty=2)

points(x, bb_data[,1]-epsi[,1], col="green", pch=".")
lines(x, bb_data[,1]-epsi[,1], col="green",lty=3)

#---------------------- Test DM  ------------------------
#---------------------- For alternative="greater", the alternative hypothesis is that method 2 is more accurate than method 1. ------------------------


forecast::dm.test(pred[,1]-epsi[,1], bb_data[,1]-epsi[,1], alternative = c("greater"), h = 1, power = 2)


