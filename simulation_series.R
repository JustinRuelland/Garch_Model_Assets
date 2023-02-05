#----------------- Simulation_series -----------------------
## Quelques fonctions utiles pour la simulation de séries financières avec le modèle GARCH(1,1)
## Il y a notamment la fonction de reconstitution des prix à partir des rendements

library(ggplot2)
library(dplyr)

###------- Reconstitution des prix à partir des rendements ----------
# On normalise le premier prix à 1
courbe_prix_arg_rendements <- function(rendements){
  n = length(rendements)
  prix = c(1) # le premier prix
  for(i in 1:n){
    prix[i+1]=prix[i]*exp(rendements[i])
  }
  m=n+1
  df = cbind(Indice = 1:m,prix)
  df = as.data.frame(df)
  
  gg_prix = ggplot(df)+geom_line(aes(x=Indice,y=prix))
  gg_prix
}


###------------------- Simulation GARCH(1,1) -----------------------
# Simulation de n carrés de rendement pour un GARCH(1,1) de paramètres theta et 
# avec un sigma carré initial

simulation_rendements <- function(n,theta, loi_eta =rnorm, sigma_init = sqrt(theta[1]/(1-theta[2]-theta[3]))){
  etas = rnorm(n)
  etas2 = etas**2
  sigmas2 = c(sigma_init**2)
  
  for(i in 2:n){
    sigmas2[i] = theta[1] + sigmas2[i-1]*(theta[2]*etas2[i-1]+theta[3])
  }
  
  epsilons2 = etas2*sigmas2
  epsilons = sign(etas)*sqrt(epsilons2)
  return(epsilons)
}

## EXEMPLE d'utilisation des deux premières fonctions
# theta = c(0.0001,0.05,0.87)
# df = simulation_rendements(1000,theta)
# courbe_prix_arg_rendements(df)


###---------------- Changement de GARCH(1,1) ------------------------
simulation_rendements_avec_changement_GARCH <- function(n,theta_1,theta_2, cut = 0.5, loi_eta = rnorm, sigma_init = sqrt(theta_1[1]/(1-theta_1[2]-theta_1[3]))){
  # Rq : je ne fais pas appel deux fois à la fonction simulation_rendements,
  # car j'ai besoin du dernier sigma2
  
  n_changement = floor(cut*n)
  
  etas = loi_eta(n)
  etas2 = etas**2
  sigmas2 = c(sigma_init**2)
  
  for(i in 2:n_changement){
    sigmas2[i] = theta_1[1] + sigmas2[i-1]*(theta_1[2]*etas2[i-1]+theta_1[3])
  }
  for(i in (n_changement+1):n){
    sigmas2[i] = theta_2[1] + sigmas2[i-1]*(theta_2[2]*etas2[i-1]+theta_2[3])
  }
  
  epsilons2 = etas2*sigmas2
  epsilons = sign(etas)*sqrt(epsilons2)
  

return(epsilons)
}

## EXEMPLE (avec theta du précis)
# theta1 = c(0.0001,0.05,0.87)
# theta2=c(0.001,0.05,0.87)
# rendements = simulation_rendements_avec_changement_GARCH(10000,theta1,theta2)
# courbe_prix_arg_rendements(rendements)


###--------------------- GARCH(2,2) --------------------------
simulation_GARCH22 <- function(n,THETA,loi_eta = rnorm, sigma_init=sqrt(THETA[1]/(1-THETA[2]-THETA[3]-THETA[4]-THETA[5]))){
  etas = rnorm(n)
  etas2 = etas**2
  
  sigmas2 = c(0.2,0.2)
  
  for (i in 3:n){
    sigmas2[i] = THETA[1] + sigmas2[i-1]*(THETA[2]*etas2[i-1]+THETA[4]) +
      sigmas2[i-2]*(THETA[3]*etas2[i-2]+THETA[5])
  }
  epsilons2 = etas2*sigmas2
  epsilons = sign(etas)*sqrt(epsilons2)
  
  return(epsilons)
}

## EXEMPLE
# THETA = c(0.0001, 0.05,0.001,0.8,0.1)
# 
# df_garch22 = simulation_GARCH22(10000,THETA)
# df_garch22_pourggplot = as.data.frame(cbind(Indice = 1:10000,df_garch22))
# 
# g = ggplot(df_garch22_pourggplot)+geom_line(aes(x = Indice,y = df_garch22))
# g
# courbe_prix_arg_rendements(df_garch22)
