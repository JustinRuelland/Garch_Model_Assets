#--------------Imports des packages--------------------
rm(list=ls())

install.packages("signal")
install.packages("tidyverse") #contient notamment ggplot2, dplyr
install.packages("lubridate")

library(signal)
library (tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate)


#---------------Import et transformation des données--------------
#setwd("C:/Users/maeld/OneDrive/Documents/GitHub/Garch_Model_Assets") #nécessaire pour Maël

source(file = "./data_preparation.R",local= TRUE)

data = read.csv("./^GDAXI.csv") #mettre le fichier de votre choix (provenant de Yahoo finance)

data <- transform_csv(data)

#Premiers graphiques - statistiques descriptives
plot_series_temp(data)
autocorrelations(data)


#--------------- Conditions de stationnarité ---------------------
source(file= "./condition_stationnarite.R",local=TRUE)

condition_stationnarite(rnorm)

runif_normalisee <-function(n){return(runif(n,-sqrt(3),sqrt(3)))}
condition_stationnarite(runif_normalisee)

rt_8_normalisee <-function(n){return(rt(n,8)/sqrt(8/(8-2)))} # Loi de Student à 8 degrés de liberté
condition_stationnarite(rt_8_normalisee)



#----------------------------QML----------------------------------------
source(file= "./QML_Mael.R",local=TRUE)

#paramètres initiaux
omega_0 <- 0.01
alpha_0 <- 0.12
beta_0 <- 0.83
theta_0 = c(omega_0,alpha_0,beta_0)
eps2_0 = 0 
sigma2_0 = omega_0/(1-alpha_0-beta_0)

#boxplots et violons
res = matrix(0,100,3)
for(i in 1:100){res[i,]=QML(simu_eps2(10**3,eps2_0,sigma2_0,theta_0))}
res = as.data.frame(res)
colnames(res) = c("omega","alpha","beta")
res$omega = res$omega-omega_0
res$alpha = res$alpha-alpha_0
res$beta = res$beta-beta_0
res_for_plot = matrix(0, ncol=2, nrow=300, byrow=FALSE)
colnames(res_for_plot) = c("value","param")
res_for_plot$value = c(res$omega,res$alpha,res$beta)
res_for_plot$param = c(rep("omega",100),rep("alpha",100),rep("beta",100))

p <- ggplot(as.data.frame(res_for_plot),aes(x=param, y=value)) + geom_boxplot()
p


p = ggplot(as.data.frame(res_for_plot), aes(x=param, y=value, fill=param)) + geom_violin(trim=FALSE)
p



#----------------normalité asymptotique----------------------
#EN COURS, sera complété avec la matrice J

#n = 10**4
#res = matrix(0,100,3)
#for(i in 1:100){res[i,]=QML(simu_eps2(n,eps2_0,sigma2_0,theta_0))$par}
#res = as.data.frame(res)
#colnames(res) = c("omega","alpha","beta")
#res$omega = sqrt(n)*(res$omega-omega_0)
#res$alpha = sqrt(n)*(res$alpha-alpha_0)
#res$beta = sqrt(n)*(res$beta-beta_0)
#plot.new() 
#hist(res$beta, breaks = 15, col = "steelblue", frame = FALSE)

#install.packages("pracma")
#library(pracma)


#s_estim= sd(res$beta)

#par(new = T)
#plot(dnorm(linspace(-5, 5, n = 100),mean=0,sd=s_estim))


#-----------------------backtest----------------------
source(file= "./prevision.R",local=TRUE)
data = read.csv("./^GDAXI.csv")
data <- transform_csv(data)
plot(data$rendement2,type='l')
vect_temp <-data$rendement2[!is.na(data$rendement2)]
eps2 = as.vector(vect_temp/sd(vect_temp)) #normalisation (REVOIR)

func_backtest(eps2,-1.96,1.96,3000) #loi normale 5%
func_backtest(eps2,sqrt(8/6)*qt(0.025,8),qt(0.975,8)*sqrt(8/6),3000) #student(8) 5%
