#-----------------Backtests--------------------

func_backtest <- function(eps,q_inf,q_sup,cut){
  
  n = length(eps)
  eps2 = eps**2
  
  #QML, on s'attend à ce que la fonction ait été importée
  theta = QML(eps2[c(1:cut)]) #estimation de theta par QML sur les données jusqu'au "cut"
  sigma2 = simu_sigma2(eps2,theta) #estimation des sigma2 à partir des epsilon2 passés
  sigma = sqrt(sigma2)
  
  #tests intervalle de confiance
  cut = cut+1
  vect_test = eps[cut:n]/sigma[cut:n]
  test_sigma = sigma[cut:n]
  outside = length(vect_test[(vect_test>q_sup)|(vect_test<q_inf)])
  inside = n - cut + 1 -outside

  #test d'adéquation du khi2 sur les Bernoulli (paramètre 0.95) (REVOIR Binomiale)
  obs = c(outside,inside) #effectifs observés
  proba = c(0.05,0.95) #probabilités théoriques
  return(chisq.test(obs,p=proba)$p.value)}
  
