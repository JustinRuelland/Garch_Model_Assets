rm(list=ls())
setwd(dir = "/Users/justinr/Documents/2A/Stat App/Github/Garch_Model_Assets/Code_Justin_R/Code_R_v2")
library(ggplot2)
library(dplyr)

source(file="./../../simulation_series.R",local=TRUE)
source(file="./../../QML_Variance.R",local=TRUE)



theta_1a <- c(0.0001, 0.12, 0.85)
theta_2a <- c(0.0001, 0.32, 0.65)
theta_3a <- c(0.0001, 0.02, 0.50)
eps1 <- simulation_rendements(2000, theta_1a)

sigma2_1 <- simu_sigma2(eps1**2, theta_1a)
sigma2_1_2a <- simu_sigma2(eps1**2, theta_2a)
sigma2_1_3a <- simu_sigma2(eps1**2, theta_3a)

dtfrm <- as.data.frame(sigma2_1, ncol=1)
dtfrm <- cbind(1:length(dtfrm[,1]), dtfrm)
dtfrm <- cbind(dtfrm,as.data.frame(sigma2_1_2a, ncol=1))
dtfrm <- cbind(dtfrm, as.data.frame(sigma2_1_3a, ncol=1))

colnames(dtfrm) <- c("Date","True_sigma2","False_sigma2_a","False_sigma2_b")

plott <- ggplot(data=dtfrm, aes(x=Date))
plott <- plott + geom_line(aes(y=True_sigma2), linetype=2, color="blue")
plott <- plott + geom_line(aes(y=False_sigma2_a),linetype=2, color="red") + ggtitle("alpha + beta constant")
#plott <- plott + guides(color="legend")
plott

plott2 <- ggplot(data=dtfrm, aes(x=Date))
plott2 <- plott2 + geom_line(aes(y=True_sigma2), linetype=2, color="blue")
plott2 <- plott2 + geom_line(aes(y=False_sigma2_b),linetype=2, color="red") + ggtitle("alpha + beta non constant")
plott2




# Erreur visualisation

erreur_moyenne <- function(serie1, serie2){
  n <- length(serie1)
  return((1/n)*sum((serie1-serie2)**2))
}

erreur_moyenne(sigma2_1, sigma2_1_2a)
erreur_moyenne(sigma2_1, sigma2_1_3a)

alpha <- c()
beta <- c()
em <- c()
truth_theta <- c(0.0001, 0.47, 0.5)
eps2 <- simulation_rendements(2000, theta_1a)**2
simu_truth <- simu_sigma2(eps2, truth_theta)

for(i in seq(0,1,by = 0.001)){
  for(j in seq(0,1,by=0.001)){
    test_theta <- c(0.0001, i, j)
    alpha <- c(alpha,i)
    beta <- c(beta,j)
    em_potential <- erreur_moyenne(simu_truth, simu_sigma2(eps2, test_theta))
    if(em_potential > 1e-05){
      em_potential <- 1e-05
    }
    em <- c(em, em_potential)
  }
  if((i*100)%%1 == 0){
    cat(i/0.01,"% \t")
  }
}

em <- (1e-05)-em
df = as.data.frame(cbind(alpha,beta,em))
p = ggplot(data = df,aes(x = alpha,y=beta,weight=em))+geom_bin2d( )+geom_point(aes(x=truth_theta[2], y=truth_theta[3]), colour="red")
p



