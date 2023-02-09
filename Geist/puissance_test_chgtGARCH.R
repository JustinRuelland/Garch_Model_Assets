###--------------- Puissance du test --------------------
# Pas le calcul de la puissance du test, 
# mais vérification pour des cas précis que le test rejette bien
rm(list=ls())

library("dplyr")
library("purrr")
library("plyr")
library("ggplot2")

source(file="./simulation_series.R")
source(file= "./QML_Variance.R")
source(file= "./prevision.R",local=TRUE)
#--------------------------- Rprof test --------------------------------

#Rprof(filename = "Rprof_carte_bleue.out", append = FALSE, interval = 0.02,
#      +       memory.profiling = FALSE, gc.profiling = FALSE,
#      +       line.profiling = FALSE, filter.callframes = FALSE,
#      +       numfiles = 100L, bufsize = 10000L)
Rprof()


#--------------------------- Changement de GARCH ------------------------------
niveau_test = 0.05
n = 1000 #nombres de jours
n_path = 25 #nombres de trajectoires simulées pour chaque carré

theta1 = c(0.0001,0.12,0.85)

N = 30 #nombre de alphas
M = 30 #nombre de betas 

a = seq(0.0,1,length.out=N)
b = seq(0,1,length.out=M)
couple = cross2(a,b)


test_for_alphabeta <- function(alpha_beta){
  rendements = simulation_rendements_avec_changement_GARCH(n,theta1,unlist(c(0.0001,c(alpha_beta))),0.8)
  p_val = func_backtest(rendements,-1.96,1.96,8*n%/%10,FALSE)$p.value #8*n%/%10 : le cut est à 80% des données
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


p = ggplot(data = df,aes(x = tab_alpha,y=tab_beta,weight=tab_res_test))+geom_bin2d( )
p

Rprof(NULL)
summaryRprof()