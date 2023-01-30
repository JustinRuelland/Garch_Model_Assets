#-----------------Backtests--------------------

func_backtest <- function(eps,q_inf,q_sup,cut,empirical){
  
  n = length(eps)
  eps2 = eps**2
  
  #QML, on s'attend à ce que la fonction ait été importée
  theta = QML(eps2[c(1:cut)]) #estimation de theta par QML sur les données jusqu'au "cut"
  sigma2 = simu_sigma2(eps2,theta) #estimation des sigma2 à partir des epsilon2 passés
  sigma = sqrt(sigma2)
  eta_quantile = eps[1:cut]/sigma[1:cut]
  
  #tests intervalle de confiance
  cut = cut+1
  eps_test = eps[cut:n]
  sigma_test = sigma[cut:n]
  
  
  
  if(empirical==FALSE){
    upper = q_sup*sigma_test
    lower = q_inf*sigma_test}
  else {upper = quantile(x=eta_quantile,prob=0.975)*sigma_test
  lower = quantile(x=eta_quantile,prob=0.025)*sigma_test}
  
  outside = length(eps_test[(eps_test>upper)|(eps_test<lower)])
  inside = n - cut + 1 -outside

  #test d'adéquation du khi2 sur les Bernoulli (paramètre 0.95)
  obs = c(outside,inside) #effectifs observés
  proba = c(0.05,0.95) #probabilités théoriques
  
  res = matrix(0,ncol=3)
  colnames(res) = c("p.value","upper.bounds","lower.bounds")
  res$p.value = chisq.test(obs,p=proba)$p.value
  res$upper.bounds = upper
  res$lower.bounds = lower

  return(res)}

#hypothèse: les etas suivent une loi normale
prevision_square <-function(eps2,cut){
  
  n = length(eps2)
  theta = QML(eps2[c(1:cut)]) #estimation de theta par QML sur les données jusqu'au "cut"
  sigma2_hat = simu_sigma2(eps2,theta)
  
  cut = cut+1
  eps2_test = eps2[cut:n]
  sigma2_test = sigma2_hat[cut:n]
  upper = qchisq(df=1,p=0.95)*sigma2_test
  outside = length(eps2_test[(eps2_test>upper)])
  inside = n - cut + 1 -outside
  
  #test d'adéquation du khi2 sur les Bernoulli (paramètre 0.95)
  obs = c(outside,inside) #effectifs observés
  proba = c(0.05,0.95) #probabilités théoriques
  
  res = matrix(0,ncol=2)
  colnames(res) = c("p.value","upper.bounds")
  res$p.value = chisq.test(obs,p=proba)$p.value
  res$upper.bounds = upper
  
  return(res)}
