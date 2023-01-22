library(dplyr)

#-----------------Backtests--------------------

func_backtest <- function(eps2,q_inf,q_sup,cut){
  
  #QML
  source(file = "C:/Users/maeld/OneDrive/Documents/GitHub/Garch_Model_Assets/QML_Mael.R",local= TRUE)
  n = length(eps2)
  theta = QML(eps2[c(1:cut)]) #estimation de theta par QML sur les données jusqu'au "cut"
  sigma2 = simu_sigma2(eps2,theta) #estimation des sigma2 à partir des epsilon2 passés
  
  #tests intervalle de confiance
  cut = cut+1
  vect_test = sqrt(eps2/sigma2)[cut:n]
  outside = length(vect_test[(vect_test>q_sup)|(vect_test<q_inf)]) #REVOIR
  inside = length(vect_test) - outside

  #test d'adéquation du khi2 sur les Bernoulli (paramètre 0.95) (REVOIR Binomiale)
  obs = c(outside,inside) #effectifs observés
  proba = c(0.05,0.95) #probabilités théoriques
  return(chisq.test(obs,p=proba)$p.value)}
  
#---------------récupération des données------------------
data = read.csv("C:/Users/maeld/OneDrive/Documents/Ensae/2A/Statapp/donnees/GOOG.csv")
source(file = "C:/Users/maeld/OneDrive/Documents/GitHub/Garch_Model_Assets/data_preparation.R",local= TRUE)
data <- transform_csv(data)
plot(data$rendement2)
eps2 = as.vector(data$rendement2/sd(data$rendement2)) #normalisation (REVOIR)

func_backtest(eps2,-1.96,1.96,200) #loi normale 5%
func_backtest(eps2,0,3.84,200) #khi2 5%
