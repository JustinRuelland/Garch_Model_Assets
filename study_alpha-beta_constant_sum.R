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
sigma2_hat_alpha_beta_constant <- function(alpha2=0.9, n=1000,seed = 8, alpha_diff_sum = 0.01, beta_diff_sum = 0.97){
  #--------------- initialisation des paramètres -----------------------
  
  # theta1 et theta2 font ici référence aux deux theta étudiés pour la deuxième partie

  cut = 0.8 #cut du backtest ET du changement de GARCH
  
  omega = 0.0001
  
  theta0 = c(omega, 0.12, 0.85)
  
  theta1 = theta0 # pas de changement de la spécification du GARCH(1,1) qui génère les données
  theta2 = c(omega, alpha2,theta1[2]+theta1[3]-alpha2) # changement, mais qui reste sur la droite alpha+beta constant
  theta_diff_sum = c(omega, alpha_diff_sum, beta_diff_sum) # changement, qui est en dehors de la droite alpha+beta constant
  
  #------------------- simulation ------------------------
  set.seed(seed)
  etas = rnorm(n)
  
  eps_1 = simulation_rendements_avec_changement_GARCH(n,theta_1 = theta0, theta_2 = theta1, cut = cut, etas = etas)
  eps_2 = simulation_rendements_avec_changement_GARCH(n,theta_1 = theta0, theta_2 = theta2, cut = cut, etas = etas)
  eps_3 = simulation_rendements_avec_changement_GARCH(n,theta_1 = theta0, theta_2 = theta_diff_sum, cut = cut, etas = etas)
  
  # Affichage des coubres de prix
  print(courbe_prix_arg_rendements(eps_1),title="Pas de changement")
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
  theta3_hat = QML(eps_3[1:n_cut]**2)
  
  sigma1_hat = simu_sigma2(eps_1**2, theta=theta1_hat)
  sigma2_hat = simu_sigma2(eps_2**2, theta=theta2_hat)
  sigma3_hat = simu_sigma2(eps_3**2, theta=theta3_hat)
  
  # Récupération des vrais sigma2
  sigma1 = eps_1**2/etas**2
  sigma2 = eps_2**2/etas**2
  sigma3 = eps_3**2/etas**2
  
  # Création du graphique final (que de n_cut à n -> test set du backtest)
  window = n_cut:n
  
  df = cbind(Date = window,sigma1_hat[window],sigma2_hat[window])
  df = as.data.frame(df)
  colnames(df) = c("Date","sigma1_hat","sigma2_hat")
  
  
  df_sigma1_hat = as.data.frame(cbind(Date = window, sigma1_hat[window]))
  df_sigma2_hat = as.data.frame(cbind(Date = window, sigma2_hat[window]))
  
  df_sigma1 = as.data.frame(cbind(Date = window, sigma1[window]))
  df_sigma2 = as.data.frame(cbind(Date = window, sigma2[window]))
  df_sigma3 = as.data.frame(cbind(Data = window, sigma3[window]))
  
  # Dataframe pour plot avec des couleurs différentes
  df_ggplot = bind_rows("Estimation sans changement" = df_sigma1_hat, "Estimation pour changement" = df_sigma2_hat, "Réel sans changement"=df_sigma1,"Réel avec changement"=df_sigma2, .id = "Courbes")
  
  print(paste0("Cas sans changement : alpha+beta = ",theta1[2]+theta1[3]))
  print(paste0("Cas avec changement (sur droite) : alpha+beta = ",theta2[2]+theta2[3]))
  print(paste0("Cas avec changement (hors droite) : alpha+beta =", theta_diff_sum[2]+theta_diff_sum[3]))
  
  p = ggplot(df_ggplot) + geom_line(aes(x=Date, y = V2, color = Courbes)) + ylab("Sigma carré") + ggtitle("Comparaison pour alpha+beta constant")#+ggtitle("Graphique des sigma carré en fonction du temps\nsans changement de GARCH et avec changement \n(pour alpha+beta constant)")
  p
  ggsave("Study_alpha2 = xx.png",width=10,height=5,path="./Graphiques_pour_Latex/Study_ab_constant_sum/")
  
  print("-------------")
  print(paste0("omega/(1-alpha-beta) = ",theta1[1]/(1-theta1[2]-theta1[3])))
  print("")
  print(paste0("Moyenne des sigma2 réels sans changement : ",mean(sigma1[n_cut:n])))
  print(paste0("Moyenne des sigma2 réels avec changement : ",mean(sigma2[n_cut:n])))
  print(paste0("Moyenne des sigma2 réels avec changement (hors droite): ",mean(sigma3[n_cut:n])))
  print("")
  print(paste0("Moyenne des sigma2 estimés sans changement : ",mean(sigma1_hat[n_cut:n])))
  print(paste0("Moyenne des sigma2 estimés avec changement : ",mean(sigma2_hat[n_cut:n])))
  print(paste0("Moyenne des sigma2 estimés avec changement (hors droite): ",mean(sigma3_hat[n_cut:n])))
  
  # Graph 2 : Sigma carré réel
  df_sigma1 = as.data.frame(cbind(Date = window, sigma1[window]))
  df_sigma2 = as.data.frame(cbind(Date = window, sigma2[window]))
  df_sigma3 = as.data.frame(cbind(Date = window, sigma3[window]))
  
  df_ggplot_reel = bind_rows("sans changement"=df_sigma1,"avec changement"=df_sigma2,"avec changement hors\n de la droite"=df_sigma3, .id = "Courbes")
  plot_reel = ggplot(df_ggplot_reel) + geom_line(aes(x = Date, y=V2, color= Courbes))+ ylab("Sigma carré") +ggtitle("Graphe des sigma2 réels")
  ggsave("Graphe_sigma2 réels.png",width=10,height=5,path="./Graphiques_pour_Latex/Study_ab_constant_sum/")
  
  print(plot_reel)
  return(p)
}


moyennes_sigma2_in_test <-function (alpha_chgt = 0.9){
  n_path = 5
  cut = 0.8
  
  omega = 0.0001
  theta1 = c(omega, 0.12, 0.85)
  theta2 = c(omega,alpha_chgt,0.97-alpha_chgt)

  res = matrix(0, nrow = 2, ncol = 3)
  rownames(res)=paste0("Changement",c("Faux","Vrai"))
  colnames(res)=paste0("n=",seq(3,5))
  
  
  for(exp in seq(3,5)){
    n = 10**exp
    for(i in 1:n_path){
      

      etas = rnorm(n) #pour avoir les mêmes etas dans les deux simulations (avec ou sans changement)
      
      eps_1 = simulation_rendements_avec_changement_GARCH(n,theta_1 = theta1, theta_2 = theta1, cut = cut, etas = etas)
      eps_2 = simulation_rendements_avec_changement_GARCH(n,theta_1 = theta1, theta_2 = theta2, cut = cut, etas = etas)
      
      n_cut = floor(cut*n)
      theta1_hat = QML(eps_1[1:n_cut]**2)
      theta2_hat = QML(eps_2[1:n_cut]**2)
      
      sigma1_hat = simu_sigma2(eps_1**2, theta=theta1_hat)
      sigma2_hat = simu_sigma2(eps_2**2, theta=theta2_hat)
      
      sigma1 = eps_1**2/etas**2
      sigma2 = eps_2**2/etas**2
      
      res[1,exp-2] =res[1,exp-2] + mean(sigma1[(0.8*n):n])
      res[2,exp-2] =res[2,exp-2] + mean(sigma2[(0.8*n):n])
    }
  }
  
  res = res/n_path
  
  return(res)
  
}

test_puissance_changment_horizon_long <- function (n=10000, alpha_chgt = 0.9){ # le but de cette fonction est de confirmé que la droite sur la "carte bleue" provient d'une vitesse de convergence faible
  loi_eta = rnorm
  
  # Paramètres GARCH
  theta1 = c(0.0001,0.12,0.85)
  beta_chgt = theta1[2]+theta1[3]-alpha_chgt
  
  # Paramètres tests
  niveau_test = 0.05
    
  test_chgt <- function(){
    rendements = simulation_rendements_avec_changement_GARCH(n,theta1,unlist(c(0.0001,alpha_chgt,beta_chgt)),0.8,etas = loi_eta(n))
    p_val = func_backtest(rendements,-1.96,1.96,0.8,TRUE)$p.value
    return(p_val)
  }
  
  n_path = 100
  sum = 0
  for(i in 1:n_path){
    if(test_chgt()<niveau_test){
      sum = sum + 1
    }
    
  }
  
  return(sum/n_path)
}
