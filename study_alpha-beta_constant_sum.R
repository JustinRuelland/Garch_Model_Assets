# In the work on the power of the test for GARCH change (puissance_test_chgtGARCH.R), 
# we found (in the case where the change in training set and test set corresponds to the date of GARCH change), 
# the existence of a straight line with slope -1 where all pairs (alpha,beta) for the second GARCH 
# were not rejected by the chi-square test.

# Hence, we focus here on the understanding of the model for two pairs of (alpha,beta) 
# such that alpha+beta is constant.
rm(list=ls())

#---------------- import -------------------
source(file = "./simulation_series.R", local=TRUE)
source(file = "./QML_Variance.R", local=TRUE)
source(file= "./prevision.R", local=TRUE)


#-------------- Studied pairs (alpha,beta) ---------------

#alpha2 must be included in [0,0,97] (0.97 = alpha1+beta1)
sigma2_hat_alpha_beta_constant <- function(alpha2=0.9, n=1000,seed = 8){
  #--------------- initialisation des paramètres -----------------------
  
  # theta1 et theta2 font ici référence aux deux theta étudiés pour la deuxième partie

  cut = 0.8 #cut du backtest ET du changement de GARCH
  
  omega = 0.0001
  
  theta0 = c(omega, 0.12, 0.85)
  
  theta1 = theta0
  theta2 = c(omega,alpha2,theta1[2]+theta1[3]-alpha2)
  
  #------------------- simulation ------------------------
  set.seed(seed)
  etas = rnorm(n)
  
  eps_1 = simulation_rendements_avec_changement_GARCH(n,theta_1 = theta0, theta_2 = theta1, cut = cut, etas = etas)
  eps_2 = simulation_rendements_avec_changement_GARCH(n,theta_1 = theta0, theta_2 = theta2, cut = cut, etas = etas)
  
  # Affichage des coubres de prix
  print(courbe_prix_arg_rendements(eps_1))
  print(courbe_prix_arg_rendements(eps_2))
  
  #--------------------- STUDY -------------------
  # Calcul des p-valeurs (pour comprendre - à ENLEVER à TERME)
  pval_1 = func_backtest(eps_1,-1.96,1.96,empirical = FALSE, cut = 0.8)$p.val
  pval_2 = func_backtest(eps_2,-1.96,1.96,empirical = FALSE, cut = 0.8)$p.val
  
  print(pval_1)
  print(pval_2)
  
  # Estimations
  n_cut = floor(cut*n)
  theta1_hat = QML(eps_1[1:n_cut]**2)
  theta2_hat = QML(eps_2[1:n_cut]**2)
  
  sigma1_hat = simu_sigma2(eps_1**2, theta=theta1_hat)
  sigma2_hat = simu_sigma2(eps_2**2, theta=theta2_hat)
  
  # Récupération des vrais sigma2
  sigma1 = eps_1**2/etas**2
  sigma2 = eps_2**2/etas**2
  
  # Création du graphique final (que de n_cut à n -> test set du backtest)
  window = n_cut:900
  
  df = cbind(Date = window,sigma1_hat[window],sigma2_hat[window])
  df = as.data.frame(df)
  colnames(df) = c("Date","sigma1_hat","sigma2_hat")
  
  
  df_sigma1_hat = as.data.frame(cbind(Date = window, sigma1_hat[window]))
  df_sigma2_hat = as.data.frame(cbind(Date = window, sigma2_hat[window]))
  
  df_sigma1 = as.data.frame(cbind(Date = window, sigma1[window]))
  df_sigma2 = as.data.frame(cbind(Date = window, sigma2[window]))
  
  # Dataframe pour plot avec des couleurs différentes
  df_ggplot = bind_rows("Estimation sans changement" = df_sigma1_hat, "Estimation pour changement" = df_sigma2_hat, "Réel sans changement"=df_sigma1,"Réel avec changement"=df_sigma2, .id = "Courbes")
  
  print(theta1[2]+theta1[3])
  print(theta2[2]+theta2[3])
  
  p = ggplot(df_ggplot) + geom_line(aes(x=Date, y = V2, color = Courbes)) + ylab("Sigma carré")#+ggtitle("Graphique des sigma carré en fonction du temps\nsans changement de GARCH et avec changement \n(pour alpha+beta constant)")
  p
  ggsave("Study_alpha2 = xx.png",width=10,height=5,path="./Graphiques_pour_Latex/Study_ab_constant_sum/")
  
  print(paste0("Moyenne sans changement : ",mean(sigma1_hat[n_cut:n])))
  print(paste0("Moyenne avec changement :",mean(sigma2_hat[n_cut:n])))

  
  return(p)
  }

sigma2_hat_alpha_beta_constant(alpha2 = 0.8, seed=31)
# Problème : la moyenne des sigma2 (les vrais déjà !!) de la série 
# avec changement est toujours environ la moitiée de la série sans changement