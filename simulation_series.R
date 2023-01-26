#----------------- Simulation_series -----------------------
## Quelques fonctions utiles pour la simulation de séries financières avec le modèle GARCH(1,1)
## Il y a notamment la fonction de reconstitution des prix à partir des rendements

library(ggplot2)

###----- Reconstitution des prix à partir des rendements ----------
# On normalise le premier prix à 1
courbe_prix_w_rendements <- function(rendements){
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


###------------------------------------------
# Simulation de n carrés de rendement pour un GARCH(1,1) de paramètres theta et 
# avec un rendement au carré initial

simulation_epsilon <- function(n,theta, loi_eta =rnorm, sigma_init = sqrt(theta[1]/(1-theta[2]-theta[3]))){
  
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

# Exempe d'utilisation des deux premières fonctions
theta = c(0.0001,0.05,0.87)
df = simulation_epsilon(1000,theta)
courbe_prix_w_rendements(df)







# A FAIRE
###-------------------------------------
simulation_epsilon_avec_changement_GARCH <- function(n,theta_1,theta_2){
  n_changement = n%/%2
  
}


simu_eps2(1000,0.06,0.5,c(0.01,0.1,0.5))
