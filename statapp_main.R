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
setwd("C:/Users/maeld/OneDrive/Documents/GitHub/Garch_Model_Assets") #nécessaire pour Maël

source(file = "./data_preparation.R",local= TRUE)

data = read.csv("./CAC40_15_19.csv") #fichier csv de Yahoo finance

data <- transform_csv(data)

#Premiers graphiques - statistiques descriptives
plot_series_temp(data)
autocorrelations(data)


#--------------- Conditions de stationnarité ---------------------
source(file= "./condition_stationnarite.R",local=TRUE)

# Définition des lois normalisées pour les eta qu'on veut tester
runif_normalisee <-function(n){return(runif(n,-sqrt(3),sqrt(3)))}
rt_8_normalisee <-function(n){return(rt(n,8)/sqrt(8/(8-2)))} # Loi de Student à 8 degrés de liberté

# Affichage et sauvergarde (dans le dossier "Graphiques_pour_Latex") du graphique qui contient les courbes de condition de stationnarité pour les 3 lois spécifiées
superposition_3graphiques_condi_statio(condition_stationnarite(rnorm),condition_stationnarite(runif_normalisee),condition_stationnarite(rt_8_normalisee))


#----------------------------QML----------------------------------------
source(file= "./QML_Variance.R",local=TRUE)

#paramètres initiaux
omega_0 <- 0.1
alpha_0 <- 0.12
beta_0 <- 0.83
theta_0 = c(omega_0,alpha_0,beta_0)
eps2_0 = 0 
sigma2_0 = omega_0/(1-alpha_0-beta_0)

#boxplots et violons
res = matrix(0,100,3)
for(i in 1:100){res[i,]=QML(simu_eps(10**3,eps2_0,sigma2_0,theta_0)**2)}
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

#QML sur le CAC40 en 2015/2016
eps2_cac = data$rendement2[1:500]
plot(eps2_cac,type='l')
QML(eps2_cac)

#tests d'invariance de alpha et beta
QML(1000*eps2_cac)


#----------------normalité asymptotique----------------------
#A reprendre, très vite fait

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


#---------------matrice de variance asymptotique-------------
omega_0 <- 0.001
alpha_0 <- 0.12
beta_0 <- 0.83
theta_0 = c(omega_0,alpha_0,beta_0)
eps2_0 = 0 
sigma2_0 = omega_0/(1-alpha_0-beta_0)

#estimation de la matrice de variance asymptotique (via le Th ergodique)
n = 10**3
eps2_sim =simu_eps(n,eps2_0,sigma2_0,theta_0)**2
var_estim = var_asymp(eps2_sim)

#estimation des coefficients diagonaux avec le RMSE
N = 100
res = matrix(0,N,3)
for(i in 1:N){res[i,]=QML(simu_eps(n,eps2_0,sigma2_0,theta_0)**2)}
res = as.data.frame(res)
colnames(res) = c("omega","alpha","beta")

#comparaison
c(sqrt(var_estim[1,1]),sqrt(var_estim[2,2]),sqrt(var_estim[3,3]))/sqrt(n)
c(sd(res$omega),sd(res$alpha),sd(res$beta))


#matrice de variance asymptotique sur le cac40
eps2_cac = data$rendement2[1:500]
var_asymp(eps2_cac)

#--------------- Puissance du test --------------------
# Pas le calcul de la puissance du test, 
# mais vérification pour des cas précis que le test rejette bien 
# source(file="./simulation_series.R")



#-----------------------backtest----------------------

source(file= "./prevision.R",local=TRUE)
source(file= "./QML_Variance.R",local=TRUE) # à importer si pas déjà fait

#a) rendements

#série simulée
theta_0 = c(10**(-4),0.12,0.83)
eps2_0 = 0 
sigma2_0 = theta_0[1]/(1-theta_0[2]-theta_0[3])
n = 10**4
eps_sim =simu_eps(n,eps2_0,sigma2_0,theta_0)
res = func_backtest(eps_sim,-1.96,1.96,6000,empirical=FALSE)
print(res$p.value)

x = c(1:4000)
plot(x, eps_sim[6001:length(eps_sim)], type = "l")
lines(x, eps_sim[6001:length(eps_sim)], col = "blue")
lines(x, res$upper.bounds, col = "red")
lines(x, res$lower.bounds, col = "red")


#cac40
eps_cac = data$rendement
res_cac = func_backtest(eps_cac,-1.96,1.96,600,empirical=FALSE) #loi normale 5%
res_cac_emp = func_backtest(eps_cac,-1.96,1.96,600,empirical=TRUE)
print(res_cac$p.value)
print(res_cac_emp$p.value)

x = c(1:421)
plot(x, eps_cac[601:length(eps_cac)], type = "l")
lines(x, eps_cac[601:length(eps_cac)], col = "blue")
lines(x, res_cac$upper.bounds, col = "red")
lines(x, res_cac$lower.bounds, col = "red")
lines(x, res_cac_emp$upper.bounds, col = "green")
lines(x, res_cac_emp$lower.bounds, col = "green")



#b) carré des rendements

#série simulée
theta_0 = c(10**(-4),0.12,0.83)
eps2_0 = 0 
sigma2_0 = theta_0[1]/(1-theta_0[2]-theta_0[3])
n = 3*10**3
eps2_sim = simu_eps(n,eps2_0,sigma2_0,theta_0)**2
res = prevision_square(eps2_sim,2000)
print(res$p.value)

x = c(1:1000)
plot(x, eps2_sim[2001:length(eps2_sim)], type = "l")
lines(x, eps2_sim[2001:length(eps2_sim)], col = "blue")
lines(x, res$upper.bounds, col = "red")  

#cac40
eps2_cac = data$rendement2
res_cac = prevision_square(eps2_cac,700) 
print(res_cac$p.value)

x = c(1:321)
plot(x, eps2_cac[701:length(eps2_cac)], type = "l")
lines(x, eps2_cac[701:length(eps2_cac)], col = "blue")
lines(x, res_cac$upper.bounds, col = "red")



#c) test de la puissance (changement de la loi des etas)

theta = c(10**(-3),0.12,0.83)
eps2_0 = 0 
sigma2_0 = theta[1]/(1-theta[2]-theta[3])
n = 10**4

#loi de laplace
rlaplace = function(n,a,b){
  u=runif(n,min=-0.5,max=0.5)
  x = a-b*sign(u)*log(1-2*abs(u))
  return(x)}

eta_lap = rlnorm(n,meanlog= 0, sdlog= 1) #loi log-normale
eta_lap = (eta_lap-mean(eta_lap))/sd(eta_lap)
print(quantile(eta_lap,probs=c(0.025,0.975)))
eta2_lap = eta_lap**2
sigmas2 = c(sigma2_0)
for(i in 2:n){sigmas2[i] = theta[1] + (theta[2]*eta2_lap[i-1] + theta[3])*sigmas2[i-1]} #sigma
eps2_lap = eta2_lap*sigmas2 
eps2_lap[0] = eps2_0
eps_lap = sign(eta_lap)*sqrt(eps2_lap)

#quantiles de la loi normale
res_lap = func_backtest(eps_lap,-1.96,1.96,6000,empirical=FALSE)
res_lap_emp = func_backtest(eps_lap,-1.96,1.96,6000,empirical=TRUE)
print(res_lap$p.value)
print(res_lap_emp$p.value)
x = c(1:4000)
plot(x, eps_lap[6001:length(eps_lap)], type = "l")
lines(x, eps_lap[6001:length(eps_lap)], col = "blue")
lines(x, res_lap$upper.bounds, col = "red")
lines(x, res_lap$lower.bounds, col = "red")
lines(x, res_lap_emp$upper.bounds, col = "green")
lines(x, res_lap_emp$lower.bounds, col = "green")
