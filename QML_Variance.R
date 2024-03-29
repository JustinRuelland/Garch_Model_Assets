#---------------------- calcul des sigma et des eps ------------------------

#calcul récurssif pour un sigma(t) précis
func_sigma2 <- function(t,sigma2_init,eps2,theta){
  res <- sigma2_init
  if(t>1){
    res<-theta[1] + theta[2]*eps2[t-1] + theta[3]*func_sigma2(t-1,sigma2_init,eps2,theta)}
  return(res)}

#simulation d'une suite de n epsilon à partir des sigma aussi simulés
simu_eps <- function(n,eps2_init,sigma2_init,theta){
  
  eta = rnorm(n,mean=0,sd=1)
  eta2 = eta**2
  
  sigmas2 = c(sigma2_init)
  for(i in 2:n){sigmas2[i] = theta[1] + (theta[2]*eta2[i-1] + theta[3])*sigmas2[i-1]} #sigma
  
  eps2 = eta2*sigmas2 
  eps2[0] = eps2_init
  eps = sign(eta)*sqrt(eps2)
  
  return(eps)}

#pareil que précedemment mais avec ue student pour les etas
simu_eps_stud <- function(n,eps2_init,sigma2_init,theta){
  
  eta = rt(n,df=8)/sqrt(8/6)
  eta2 = eta**2
  
  sigmas2 = c(sigma2_init)
  for(i in 2:n){sigmas2[i] = theta[1] + (theta[2]*eta2[i-1] + theta[3])*sigmas2[i-1]} #sigma
  
  eps2 = eta2*sigmas2 
  eps2[0] = eps2_init
  eps = sign(eta)*sqrt(eps2)
  
  return(eps)}

#simulation des sigmas**2 à partir des epsilon**2 connus
simu_sigma2 <- function(eps2,theta){
  n = length(eps2)
  
  sigmas2 = double(n) #création d'un vecteur de taille qui ne va pas varier par la suite
  sigmas2[1] = c(mean(eps2))
  
  for(i in 2:n){sigmas2[i]=theta[1]+theta[2]*eps2[i-1]+theta[3]*sigmas2[i-1]}
  
  return(sigmas2)
}


#---------------------- Methode du quasi-maximum de vraisemblance ----------

#on part du principe qu'on a observe une série financière
#dont les rendements au carré sont eps2
QML <-function(eps2){
  n = length(eps2) #longueur de la série financière
  
  #fonction de vraisemblance à optimiser
  f_opt <- function(theta_opt){
    sigmas2_QML = simu_sigma2(eps2,theta_opt) #simulation des sigmas à partir des eps2
    #on retire les 20 premières données (négligeables cf notes)
    return(sum(log(sigmas2_QML[25:n])+(eps2[25:n]/sigmas2_QML[25:n]))) } #log vraisemblance
  
  theta_init = c(0.0001,0.18,0.8) #valeur initiale
  
  #contraintes
  ui <- cbind(c(1,-1,0,0,0,0),c(0,0,1,-1,0,0),c(0,0,0,0,1,-1))
  ci <- c(10**(-9), -1, 0, -3, 0, -0.99)
  return(constrOptim(theta=theta_init,f = f_opt,ci=ci,ui=ui,gr=NULL)$par)} #opti sous contraintes linéaires


#-------------estimation de la matrice de variance asymptotique-----------------------------

var_asymp <- function(eps2){
  
  n = length(eps2)
  
  #matrice J
  
  #calcul itéré du gradient
  iter_grad <- function(grad,theta,eps2_t,sigma2_t){
    new_grad = c()
    new_grad[1] = 1 + theta[3]*grad[1]
    new_grad[2] = eps2_t + theta[3]*grad[2]
    new_grad[3] = sigma2_t + theta[3]*grad[3]
    return(new_grad)}
  
  #estimation de theta par QML
  theta_estim = QML(eps2)
  
  #gradient initial
  grad = c(0,0,0)
  
  #estimation des sigma2 (sigma2_hat)
  sigma2_estim = simu_sigma2(eps2,theta_estim)
  J = (1/sigma2_estim[1]**2)*(grad%*%t(grad))/n
  for(i in 2:n){
    grad = iter_grad(grad,theta_estim,eps2[i-1],sigma2_estim[i-1])
    J = J + (1/sigma2_estim[i]**2)*(grad%*%t(grad))/n  }
  
  #coefficient K
  eta4_estim = (eps2/sigma2_estim)**2
  K = mean(eta4_estim)
  
  #formule de la variance asymptotique
  var_asymp = (K-1)*solve(J)
  
  return(var_asymp)}

