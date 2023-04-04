#-----------------Backtests--------------------

func_backtest <- function(eps,q_inf,q_sup,cut,empirical){
  
  n = length(eps)
  n_cut = floor(cut*n)
  eps2 = eps**2
  
  #QML, on s'attend à ce que la fonction ait été importée
  theta = QML(eps2[c(1:n_cut)]) #estimation de theta par QML sur les données jusqu'au "cut"
  sigma2 = simu_sigma2(eps2,theta) #estimation des sigma2 à partir des epsilon2 passés
  sigma = sqrt(sigma2)
  eta_quantile = eps[1:n_cut]/sigma[1:n_cut]
  
  #tests intervalle de confiance
  n_cut = n_cut+1
  eps_test = eps[n_cut:n]
  sigma_test = sigma[n_cut:n]
  
  
  if(empirical==FALSE){
    upper = q_sup*sigma_test
    lower = q_inf*sigma_test}
  else{
    upper = quantile(x=eta_quantile,prob=0.975)*sigma_test
    lower = quantile(x=eta_quantile,prob=0.025)*sigma_test}
  
  outside = length(eps_test[(eps_test>upper)|(eps_test<lower)])
  inside = n - n_cut + 1 -outside

  #test d'adéquation du khi2 sur les Bernoulli (paramètre 0.95)
  obs = c(outside,inside) #effectifs observés
  proba = c(0.05,0.95) #probabilités théoriques
  
  res = matrix(0,ncol=3)
  colnames(res) = c("p.value","upper.bounds","lower.bounds")
  res$p.value = chisq.test(obs,p=proba)$p.value
  res$upper.bounds = upper
  res$lower.bounds = lower

  return(res)}

#func_backtest <- function(eps,q_inf,q_sup,cut,empirical){
  
  n = length(eps)
  n_cut = floor(cut*n)
  eps2 = eps**2
  
  #QML, on s'attend à ce que la fonction ait été importée
  theta = QML(eps2[c(1:n_cut)]) #estimation de theta par QML sur les données jusqu'au "cut"
  sigma2 = simu_sigma2(eps2,theta) #estimation des sigma2 à partir des epsilon2 passés
  sigma = sqrt(sigma2)
  eta_quantile = eps[1:n_cut]/sigma[1:n_cut]
  
  #tests intervalle de confiance
  n_cut = n_cut+1
  eps_test = eps[n_cut:n]
  sigma_test = sigma[n_cut:n]
  
  
  if(empirical==FALSE){
    upper = q_sup*sigma_test
    lower = q_inf*sigma_test}
  else{
    upper = quantile(x=eta_quantile,prob=0.975)*sigma_test
    lower = quantile(x=eta_quantile,prob=0.025)*sigma_test}
  
  outside = length(eps_test[(eps_test>upper)|(eps_test<lower)])
  inside = n - n_cut + 1 -outside
  
  #test d'adéquation du khi2 sur les Bernoulli (paramètre 0.95)
  obs = c(outside,inside) #effectifs observés
  proba = c(0.05,0.95) #probabilités théoriques
  
  res = matrix(0,ncol=3)
  colnames(res) = c("p.value","upper.bounds","lower.bounds")
  res$p.value = chisq.test(obs,p=proba)$p.value
  res$upper.bounds = upper
  res$lower.bounds = lower
  
  return(res)
}

#hypothèse: les etas suivent une loi normale
backtest_square_modified <-function(eps2,q_inf, q_sup, cut, empirical){
  
  n = length(eps2)
  n_cut = floor(cut*n)
  theta = QML(eps2[c(1:n_cut)]) #estimation de theta par QML sur les données jusqu'au "cut"
  sigma2_hat = simu_sigma2(eps2,theta)
  
  n_cut = n_cut+1
  eps2_test = eps2[n_cut:n]
  sigma2_test = sigma2_hat[n_cut:n]
  
  if(empirical == F){
    upper <- qchisq(0.975,df=1)*sigma2_test
    lower <- qchisq(0.025, df=1)*sigma2_test
  }else{
    etas_square <- eps2[1:n_cut]/sigma2_hat[1:n_cut]
    upper <- quantile(etas_square, prob = 0.975)*sigma2_test
    lower <- quantile(etas_square, prob = 0.025)*sigma2_test
  }
  
  outside = length(eps2_test[(eps2_test>upper)|(eps2_test<lower)])
  inside = n - n_cut + 1 -outside
  
  
  #test d'adéquation du khi2 sur les Bernoulli (paramètre 0.95)
  obs = c(outside,inside) #effectifs observés
  proba = c(0.05,0.95) #probabilités théoriques
  
  res = matrix(0,ncol=2)
  colnames(res) = c("p.value","upper.bounds")
  res$p.value = chisq.test(obs,p=proba)$p.value
  res$upper.bounds = upper
  
  return(res)
}


backtest_square <-function(eps2,cut){
  
  n = length(eps2)
  n_cut = floor(cut*n)
  theta = QML(eps2[c(1:n_cut)]) #estimation de theta par QML sur les données jusqu'au "cut"
  sigma2_hat = simu_sigma2(eps2,theta)
  
  n_cut = n_cut+1
  eps2_test = eps2[n_cut:n]
  sigma2_test = sigma2_hat[n_cut:n]
  upper = qchisq(df=1,p=0.95)*sigma2_test
  outside = length(eps2_test[(eps2_test>upper)])
  inside = n - n_cut + 1 -outside
  
  #test d'adéquation du khi2 sur les Bernoulli (paramètre 0.95)
  obs = c(outside,inside) #effectifs observés
  proba = c(0.05,0.95) #probabilités théoriques
  
  res = matrix(0,ncol=2)
  colnames(res) = c("p.value","upper.bounds")
  res$p.value = chisq.test(obs,p=proba)$p.value
  res$upper.bounds = upper
  
  return(res)}

#------------horizon 2, estimation par noyau--------------------
simu_rejet <- function(ker){
  
  func <-function(x){return((1/sqrt(2*pi))*exp(-0.5*x**2))}
  M  = max(ker$y/func(ker$x))
  
  #pour rentrer dans la boucle
  U = 1
  V = 1
  Z = 0
  ind = 1
  cpt = 0
  
  while(U*M*func(V)>Z*ker$y[ind]){
    cpt = cpt +1
    Z = 1
    U = runif(1,0,1)
    V = rnorm(1,0,1)
    if(cpt%%2==0){ind = c(which(ker$x>V))[1]}
    else{ind = c(which(ker$x<V))[length(c(which(ker$x<V)))]}
    }
  
  return(V)}

#ker =  density(rnorm(10**4,0,1))
#c = c()
#for(i in 1:10**3){c[i]=simu_rejet(ker)}
#hist(c,breaks=80,prob=TRUE)
#mean(c)


pred_h2_kernel <- function(eps,cut,nb_sim){
  
  n = length(eps)
  n_cut = floor(cut*n)
  
  theta_hat = QML(eps[1:n_cut]**2)
  w  = theta_hat[1]
  a = theta_hat[2]
  b = theta_hat[3]
  
  sigma2 = simu_sigma2(eps**2, theta_hat)
  sigma = sqrt(sigma2)
  eta = eps[1:n_cut]/sigma[1:n_cut]
  
  ker = density(eta)

  
  n_cut = n_cut + 1
  bounds_inf = c()
  bounds_sup = c()
  
  for(i in n_cut:n){
    
    pred_h2 = c() 
    
    for(j in 1:nb_sim){
      eta2_h1 = simu_rejet(ker)
      eta2_h2 = simu_rejet(ker)
      pred_h2[j] = eta2_h2*(w + (a*eta2_h2+b)*(w + a*eps[i-2]**2 + b*sigma2[i-2]))
      
    bounds_inf[i] = quantile(x=pred_h2, probs=0.025)
    bounds_sup[i] = quantile(x=pred_h2, probs=0.975) }   }
  
  return(c(bounds_inf,bounds_sup))}



#setwd("C:/Users/maeld/OneDrive/Documents/GitHub/Garch_Model_Assets") #nécessaire pour Maël
#source(file= "./QML_Variance.R",local=TRUE)


#omega_0 <- 0.001
#alpha_0 <- 0.12
#beta_0 <- 0.83
#theta_0 = c(omega_0,alpha_0,beta_0)
#eps2_0 = 0 
#sigma2_0 = omega_0/(1-alpha_0-beta_0)
#n = 10**4
#eps_sim =simu_eps(n,eps2_0,sigma2_0,theta_0)
#pred_h2_kernel(eps_sim,8000,50)



