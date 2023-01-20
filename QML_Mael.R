library(ggplot2)
#---------------------paramètres initiaux------------
omega_0 <- 0.01
alpha_0 <- 0.12
beta_0 <- 0.82
theta_0 = c(omega_0,alpha_0,beta_0)
eps2_0 = 0 
sigma2_0 = omega_0/(1-alpha_0-beta_0)

#---------------------- calcul des sigma2 et des eps2 ------------------------

#calcul récurssif pour un sigma(t) précis (FACULTATIF)
func_sigma2 <- function(t,sigma2_init,eta2,theta){
  res <- sigma2_init
  if(t>1){
    res<-theta[1] + (theta[2]*eta2[t-1] + theta[3])*func_sigma2(t-1,sigma2_init,theta)}
  return(res)}

#calcul d'une suite de n epsilon(t)**2 (à partir des sigmas(t))
simu_eps2 <- function(n,eps2_init,sigma2_init,theta){
  eta2 = rnorm(n,mean=0,sd=1)**2
  sigmas2 = c(sigma2_init)
  for(i in 2:n){sigmas2[i] = theta[1] + (theta[2]*eta2[i-1] + theta[3])*sigmas2[i-1]}
  eps2 = eta2*sigmas2
  eps2[0] = eps2_init
  return(eps2)}

#---------------------- Methode du quasi-maximum de vraisemblance ----------
QML <-function(eps2){
  n = length(eps2)
  
  f_opt <- function(theta_opt){
    sigmas2_QML = c(theta_opt[1]/(1-theta_opt[2]-theta_opt[3]))
    for(i in 2:n){sigmas2_QML[i]=theta_opt[1]+theta_opt[2]*eps2[i-1]+theta_opt[3]*sigmas2_QML[i-1]}
    #on retire les 20 premières données (négligeables cf notes)
    return(sum(log(sigmas2_QML[100:n])+(eps2[100:n]/sigmas2_QML[100:n]))) } 
  
  theta_init = c(0.06,0.2,0.6)
  ui <- cbind(c(1,-1,0,0,0,0),c(0,0,1,-1,0,0),c(0,0,0,0,1,-1))
  ci <- c(0.001, -1, 0, -3, 0, -0.99)
  return(constrOptim(theta=theta_init,f = f_opt,ci=ci,ui=ui,gr=NULL))}

res = matrix(0,100,3)
for(i in 1:100){res[i,]=QML(simu_eps2(10**3,eps2_0,sigma2_0,theta_0))$par}
res = as.data.frame(res)
colnames(res) = c("omega","alpha","beta")
res$omega = res$omega-omega_0
res$alpha = res$alpha-alpha_0
res$beta = res$beta-beta_0
boxplot(res$omega)
boxplot(res$alpha)
boxplot(res$beta)
