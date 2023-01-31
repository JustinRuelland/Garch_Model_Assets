###--------------- Puissance du test --------------------
# Pas le calcul de la puissance du test, 
# mais vérification pour des cas précis que le test rejette bien

library("dplyr")
library("purrr")
library("plyr")
library("ggplot2")

setwd("C:/Users/louis/Garch_Model_Assets")
source(file="./simulation_series.R")
source(file= "./QML_Variance.R")
source(file= "./prevision.R",local=TRUE)

#--------------------------- Changement de GARCH --------------------------------
n = 1000 #nombres de jours
theta1 = c(0.0001,0.12,0.85)

N = 30 #nombre de alphas
M = 30 #nombre de betas 

a = seq(0.0,0.3,length.out=N)
b = seq(0.6,1,length.out=M)
couple = cross2(a,b)


f_epaule <- function(alpha_beta){
  rendements = simulation_rendements_avec_changement_GARCH(n,theta1,unlist(c(0.0001,c(alpha_beta))))
  p_val = func_backtest(rendements,-1.96,1.96,8*n%/%10,FALSE)$p.value
  return(p_val)
}

tab_alpha = c()
tab_beta = c()

long = N*M
for(i in 1:long){
  tab_alpha[i] = c(couple[[i]][[1]])
  tab_beta[i] = c(couple[[i]][[2]])
}


tab_pval = laply(couple,f_epaule)

### Echelonnage des p-valeurs
echelon_pval<-function(pval){
  if(pval<0.05){
    return(0)
  }
  else if((0.05<=pval)&(pval<0.1)){
    return(0.05)
  }
  else{
    return(1)
  }
}
# A appliquer si on veut connaître que des tranches de p-valeur
#tab_pval = laply(tab_pval,echelon_pval)


### Affichage
# Transformation données
df = as.data.frame(cbind(tab_alpha,tab_beta,tab_pval))


p = ggplot(data = df,aes(x = tab_alpha,y=tab_beta,weight=tab_pval))+geom_bin2d( )
p
