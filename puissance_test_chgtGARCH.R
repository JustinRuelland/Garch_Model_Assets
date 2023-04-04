###--------------- Puissance du test --------------------
# Pas le calcul de la puissance du test, 
# mais vérification pour des cas précis que le test rejette bien
rm(list=ls())

library("dplyr")
library("purrr")
library("plyr")
library("ggplot2")


source(file="./simulation_series.R",local=TRUE)
source(file= "./QML_Variance.R",local=TRUE)
source(file= "./prevision.R",local=TRUE)
#--------------------------- Rprof test --------------------------------
# Rprof()

#--------------------------- Changement de GARCH ------------------------------
puissance_test_chgtGARCH<-function(cut_chgt=0.4,n_path=10,n=1000,loi_eta=rnorm, backtest = func_backtest){
  #n : nombres de jours
  #n_path : nombres de trajectoires simulées pour chaque carré
  
  niveau_test = 0.05
  theta1 = c(0.0001,0.12,0.55)
  
  N = 30 #nombre de alphas
  M = 30 #nombre de betas 
  
  a = seq(0,1,length.out=N)
  b = seq(0,1,length.out=M)
  couple = cross2(a,b)
  
  
  test_for_alphabeta <- function(alpha_beta){
    rendements = simulation_rendements_avec_changement_GARCH(n,theta1,unlist(c(0.0001,c(alpha_beta))),cut_chgt,etas = loi_eta(n))
    p_val = backtest(rendements,-1.96,1.96,0.8,TRUE)$p.value #0.8 : le cut est à 80% des données
    return(p_val)
  }
  
  tab_alpha = c()
  tab_beta = c()
  
  long = N*M
  for(i in 1:long){
    tab_alpha[i] = c(couple[[i]][[1]])
    tab_beta[i] = c(couple[[i]][[2]])
  }
  
  
  tab_pval = laply(couple,test_for_alphabeta)
  
  tab_res_test = rep(0,N*M)
  for(j in 1:n_path){
    print("Numéro de trajectoire :")
    print(j)
    tab_pval = ###--------------- Puissance du test --------------------
# Pas le calcul de la puissance du test, 
# mais vérification pour des cas précis que le test rejette bien
rm(list=ls())

library("dplyr")
library("purrr")
library("plyr")
library("ggplot2")


source(file="./simulation_series.R",local=TRUE)
source(file= "./QML_Variance.R",local=TRUE)
source(file= "./prevision.R",local=TRUE)
#--------------------------- Rprof test --------------------------------
# Rprof()

#--------------------------- Changement de GARCH ------------------------------
puissance_test_chgtGARCH<-function(cut_chgt=0.4,n_path=10,n=1000,loi_eta=rnorm, backtest = func_backtest){
  #n : nombres de jours
  #n_path : nombres de trajectoires simulées pour chaque carré
  
  niveau_test = 0.05
  theta1 = c(0.0001,0.12,0.55)
  
  N = 30 #nombre de alphas
  M = 30 #nombre de betas 
  
  a = seq(0,1,length.out=N)
  b = seq(0,1,length.out=M)
  couple = cross2(a,b)
  
  
  test_for_alphabeta <- function(alpha_beta){
    rendements = simulation_rendements_avec_changement_GARCH(n,theta1,unlist(c(0.0001,c(alpha_beta))),cut_chgt,etas = loi_eta(n))
    p_val = backtest(rendements,-1.96,1.96,0.8,TRUE)$p.value #0.8 : le cut est à 80% des données
    return(p_val)
  }
  
  tab_alpha = c()
  tab_beta = c()
  
  long = N*M
  for(i in 1:long){
    tab_alpha[i] = c(couple[[i]][[1]])
    tab_beta[i] = c(couple[[i]][[2]])
  }
  
  
  tab_pval = laply(couple,test_for_alphabeta)
  
  tab_res_test = rep(0,N*M)
  for(j in 1:n_path){
    print("Numéro de trajectoire :")
    print(j)
    tab_pval = laply(couple,test_for_alphabeta)
    
    # Syntaxe pour associer la valeur 1 si la pvaleur est supérieure au niveau, 0 sinon
    tab_pval[tab_pval>niveau_test] = 1
    tab_pval[tab_pval<=niveau_test] = 0
    
    tab_res_test = tab_res_test+tab_pval
  }
  
  tab_res_test = tab_res_test /n_path # pour obtenir moyenne
  
  
  ### Affichage
  # Transformation données
  df = as.data.frame(cbind(tab_alpha,tab_beta,tab_res_test))
  
  
  p = ggplot(data = df,aes(x = tab_alpha,y=tab_beta,weight=tab_res_test))+geom_bin2d( )+geom_point(aes(x=theta1[2], y=theta1[3]), colour="red")
  print(p)
  ggsave("cb - n_path =XX - n=XXXX - cut=XX.png", width = 6, height = 5,path="./Graphiques_pour_Latex/Carte_bleue_(puissance_chgt_GARCH)/")
  return(p)
}

### Exemple
# source(file = "./white_noise_laws.R",local= TRUE)

#source(file = "./white_noise_laws.R",local=TRUE)
puissance_test_chgtGARCH(0.8,5,1000,loi_eta = rnorm, backtest = backtest_square_modified)

# Rprof(NULL)
# summaryRprof()
laply(couple,test_for_alphabeta)
    
    # Syntaxe pour associer la valeur 1 si la pvaleur est supérieure au niveau, 0 sinon
    tab_pval[tab_pval>niveau_test] = 1
    tab_pval[tab_pval<=niveau_test] = 0
    
    tab_res_test = tab_res_test+tab_pval
  }
  
  tab_res_test = tab_res_test /n_path # pour obtenir moyenne
  
  
  ### Affichage
  # Transformation données
  df = as.data.frame(cbind(tab_alpha,tab_beta,tab_res_test))
  
  
  p = ggplot(data = df,aes(x = tab_alpha,y=tab_beta,weight=tab_res_test))+geom_bin2d( )+geom_point(aes(x=theta1[2], y=theta1[3]), colour="red")
  print(p)
  ggsave("cb - n_path =XX - n=XXXX - cut=XX.png", width = 6, height = 5,path="./Graphiques_pour_Latex/Carte_bleue_(puissance_chgt_GARCH)/")
  return(p)
}

### Exemple
# source(file = "./white_noise_laws.R",local= TRUE)

#source(file = "./white_noise_laws.R",local=TRUE)
puissance_test_chgtGARCH(0.8,5,1000,loi_eta = rnorm, backtest = backtest_square_modified)

# Rprof(NULL)
# summaryRprof()
