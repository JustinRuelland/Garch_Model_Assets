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
  eps_test = eps[cut:n]
  sigma_test = sigma[cut:n]
  
  upper = q_sup*sigma_test
  lower = q_inf*sigma_test
  outside = length(eps_test[(eps_test>upper)|(eps_test<lower)])
  inside = n - cut + 1 -outside

  #test d'adéquation du khi2 sur les Bernoulli (paramètre 0.95) (REVOIR Binomiale)
  obs = c(outside,inside) #effectifs observés
  proba = c(0.05,0.95) #probabilités théoriques
  
  p = chisq.test(obs,p=proba)$p.value
  res = matrix(0,ncol=3)
  colnames(res) = c("p.value","upper.bounds","lower.bounds")
  res$p.value = p
  res$upper.bounds = upper
  res$lower.bounds = lower

  return(res)}