#--------------Imports des packages--------------------
rm(list=ls())

install.packages("signal")
install.packages("tidyverse") #contient notamment ggplot2, dplyr
install.packages("lubridate")
install.packages("forecast")

library(signal)
library (tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate)
library(zoo)
library(forecast)
library(pracma)

# DEFINITION DU WORKING DIRECTORY : dans le menu "Session", "Set working directory" 
# et choisir le dossier "\Garch_Model_Assets"

#---------------Import et transformation des données--------------
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

#----------------normalité asymptotique----------------------

n = 10**3
res = matrix(0,n,3)
V = matrix(0,3,3)
for(i in 1:n){
  eps2 = simu_eps(10**3,eps2_0,sigma2_0,theta_0)**2
  res[i,]=QML(eps2)}
res = as.data.frame(res)
colnames(res) = c("omega","alpha","beta")
res$omega = sqrt(n)*(res$omega-omega_0)
res$alpha = sqrt(n)*(res$alpha-alpha_0)
res$beta = sqrt(n)*(res$beta-beta_0)

for(i in 1:100){V = V + var_asymp(simu_eps(10**3,eps2_0,sigma2_0,theta_0)**2)/100}

a = min(res$alpha)
b = max(res$alpha)
xx = seq(from = a, to = b, by = (b-a)/50)
sd_hat = sqrt(V[2,2])
y = c()
for(i in 1:length(xx)){y[i]=dnorm(xx[i],mean=0,sd=sd_hat)}

y_hat = hist(res$alpha, breaks = xx, col = "steelblue")$density
plot(xx[2:51],y_hat,col="red",type='o')
lines(xx,y)

ks.test(res$alpha,"pnorm",mean=0, sd=sd_hat)

#--------------- Puissance du test --------------------
# Pas le calcul de la puissance du test, 
# mais vérification pour des cas précis que le test rejette bien 
# source(file="./simulation_series.R")



#-----------------------backtest----------------------

source(file= "./prevision.R",local=TRUE)
source(file= "./QML_Variance.R",local=TRUE) # à importer si pas déjà fait

#a) rendements

theta_0 = c(10**(-4),0.12,0.83)
eps2_0 = 0 
sigma2_0 = theta_0[1]/(1-theta_0[2]-theta_0[3])
n = 10**4

#série simulée, graphique avec les IC
eps_sim =simu_eps(n,eps2_0,sigma2_0,theta_0)
res = func_backtest(eps_sim,-1.96,1.96,0.8,empirical=FALSE)

#paramètres pour les plots
n_cut  = floor(0.8*length(eps_sim))
l1 = n-n_cut
l2 = n_cut+1
xx = c(1:l1)

res_for_plot = as.data.frame(res)
res_for_plot$abs= xx
res_for_plot$eps_sim = eps_sim[l2:length(eps_sim)]

ggplot(data = res_for_plot, mapping = aes(x=abs,y=eps_sim)) + geom_line(color='blue') + 
  geom_line(data = res_for_plot, aes(x=abs, y=upper.bounds), color='red') + 
            geom_line(data = res_for_plot, aes(x=abs, y=lower.bounds), color='red')

#estimation, un peu long à executer
distrib = c()
for(i in 1:1000){
  eps_sim = simu_eps(10**3,eps2_0,sigma2_0,theta_0)
  res = func_backtest(eps_sim,-1.96,1.96,0.8,empirical=FALSE)
  distrib[i] = res$p.value}

print(length(distrib[distrib<0.05])/length(distrib))



#cac40
eps_cac = data$rendement
res_cac = func_backtest(eps_cac,-1.96,1.96,0.8,empirical=FALSE) #loi normale 5%
res_cac_emp = func_backtest(eps_cac,-1.96,1.96,0.8,empirical=TRUE)
print(res_cac$p.value)
print(res_cac_emp$p.value)

n = length(eps_cac)
n_cut  = floor(0.8*n)
l1 = n-n_cut
l2 = n_cut+1
xx = c(1:l1)

res_for_plot = as.data.frame(res_cac)
res_for_plot$abs= xx
res_for_plot$eps_cac= eps_cac[l2:length(eps_cac)]
res_for_plot$upper.bounds_emp= res_cac_emp$upper.bounds
res_for_plot$lower.bounds_emp= res_cac_emp$lower.bounds

ggplot(data = res_for_plot, mapping = aes(x=abs,y=eps_cac)) + geom_line(color='blue') + 
  geom_line(data = res_for_plot, aes(x=abs, y=upper.bounds), color='red') + 
  geom_line(data = res_for_plot, aes(x=abs, y=lower.bounds), color='red') + 
  geom_line(data = res_for_plot, aes(x=abs, y=upper.bounds_emp), color='green') + 
  geom_line(data = res_for_plot, aes(x=abs, y=lower.bounds_emp), color='green')  



#b) carré des rendements

#série simulée
theta_0 = c(10**(-4),0.12,0.83)
eps2_0 = 0 
sigma2_0 = theta_0[1]/(1-theta_0[2]-theta_0[3])

eps2_sim = simu_eps(2000,eps2_0,sigma2_0,theta_0)**2
res = backtest_square(eps2_sim,0.75)

n = length(eps2_sim)
n_cut  = floor(0.75*n)
l1 = n-n_cut
l2 = n_cut+1
xx = c(1:l1)
res_for_plot = as.data.frame(res)
res_for_plot$eps2 = eps2_sim[l2:length(eps2_sim)]
res_for_plot$abs= xx
ggplot(data = res_for_plot, mapping = aes(x=abs,y=eps2)) + geom_line(color='blue') + 
  geom_line(data = res_for_plot, aes(x=abs, y=upper.bounds), color='red')

#estimation
distrib = c()
for(i in 1:1000){
  eps2_sim = simu_eps(10**3,eps2_0,sigma2_0,theta_0)**2
  res = backtest_square(eps2_sim,800)
  distrib[i] = res$p.value}

print(length(distrib[distrib<0.05])/length(distrib))

#cac40
eps2_cac = data$rendement2
res_cac = backtest_square(eps2_cac,0.75) 
print(res_cac$p.value)

n = length(eps2_cac)
n_cut  = floor(0.75*n)
l1 = n-n_cut
l2 = n_cut+1
xx = c(1:l1)
res_for_plot = as.data.frame(res_cac)
res_for_plot$eps2 = eps2_cac[l2:length(eps2_cac)]
res_for_plot$abs= xx
ggplot(data = res_for_plot, mapping = aes(x=abs,y=eps2)) + geom_line(color='blue') + 
  geom_line(data = res_for_plot, aes(x=abs, y=upper.bounds), color='red')



#c) test de la puissance 

# i. changement de la loi des etas

theta = c(10**(-3),0.12,0.83)
eps2_0 = 0 
sigma2_0 = theta[1]/(1-theta[2]-theta[3])
n = 10**3

cpt = 0
cpt_emp = 0
for(i in 1:100){
  eta_lap = runif(n,min= -1, max= 1) #loi uniforme
  eta_lap = (eta_lap-mean(eta_lap))/sd(eta_lap)
  eta2_lap = eta_lap**2
  sigmas2 = c(sigma2_0)
  for(i in 2:n){sigmas2[i] = theta[1] + (theta[2]*eta2_lap[i-1] + theta[3])*sigmas2[i-1]} #sigma
  eps2_lap = eta2_lap*sigmas2 
  eps2_lap[0] = eps2_0
  eps_lap = sign(eta_lap)*sqrt(eps2_lap)

  #quantiles de la loi normale
  res_lap = func_backtest(eps_lap,-1.96,1.96,0.8,empirical=FALSE)
  res_lap_emp = func_backtest(eps_lap,-1.96,1.96,0.8,empirical=TRUE)
  if(res_lap$p.value<0.05){cpt=cpt+1}
  if(res_lap_emp$p.value<0.05){cpt_emp=cpt_emp+1} }
  
print(c(cpt,cpt_emp)/100)

#plot de la dernière simulation pour illustrer
n = length(eps_lap)
n_cut  = floor(0.8*n)
l1 = n-n_cut
l2 = n_cut+1
xx = c(1:l1)
res_for_plot = as.data.frame(res_lap)
res_for_plot$eps2 = eps_lap[l2:length(eps_lap)]
res_for_plot$abs= xx
res_for_plot$upper.bounds_emp = res_lap_emp$upper.bounds
res_for_plot$lower.bounds_emp = res_lap_emp$lower.bounds
ggplot(data = res_for_plot, mapping = aes(x=abs,y=eps2)) + geom_line(color='blue') + 
  geom_line(data = res_for_plot, aes(x=abs, y=upper.bounds), color='red')+
  geom_line(data = res_for_plot, aes(x=abs, y=lower.bounds), color='red')+
  geom_line(data = res_for_plot, aes(x=abs, y=upper.bounds_emp), color='green')+
  geom_line(data = res_for_plot, aes(x=abs, y=lower.bounds_emp), color='green')


# ii. Changement des paramètres de GARCH(1,1)

source(file="./puissance_test_chgtGARCH.R")

puissance_test_chgtGARCH(0.4,10,1000) #Affiche la "carte bleue" - cela prend 5 à 10min avec ces paramètres


#----------------------Test de Mariano--------------------------

source(file="./comparison_model.R")
source(file= "./QML_Variance.R",local=TRUE)#si pas déjà importé
 

eps2 = data$rendement2
n = length(eps2)
n_cut = floor(0.8*n)

yy_garch = pred_h1_garch(eps2,0.8)
yy_roll = rolling_av(eps2,0.8,4,liss_exp=TRUE)#window = 4
e = n-n_cut
s = n_cut+1
res_for_plot = data.frame(abs = c(1:e), garch= yy_garch, roll= yy_roll, eps2 = eps2[s:n])
ggplot(data = res_for_plot, mapping = aes(x=abs,y=eps2)) + geom_line(color='blue') + 
  geom_line(data = res_for_plot, aes(x=abs, y=garch), color='red') + 
  geom_line(data = res_for_plot, aes(x=abs, y=roll), color='green')

test_mariano(yy_garch,yy_roll,eps2[s:n],hor=1) #test de mariano sur le cac40


#100 tests de Mariano sur des données simulées
omega_0 <- 10**(-4)
alpha_0 <- 0.12
beta_0 <- 0.83
theta_0 = c(omega_0,alpha_0,beta_0)
eps2_0 = 0 
sigma2_0 = omega_0/(1-alpha_0-beta_0)

cpt = 0
for(i in 1:100){
  eps2 = simu_eps(10**3,eps2_0,sigma2_0,theta_0)**2
  s = floor(0.8*length(eps2))+1
  t = test_mariano(pred_h1_garch(eps2,0.8),rolling_av(eps2,0.8,4,liss_exp=TRUE),eps2[s:length(eps2)],hor=1)
  if(t$statistic<0 & t$p.value<0.001){cpt = cpt+1}}

print(cpt) #nb de fois où Garch "l'emporte" (avec un niveau à 0.1%) sur 100 simulations

