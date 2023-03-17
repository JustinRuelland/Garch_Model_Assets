rm(list=ls())
setwd(dir = "/Users/justinr/Documents/2A/Stat App/Github/Garch_Model_Assets/Code_Justin_R/Code_R_v2")
library(ggplot2)
library(dplyr)

source(file="./../../simulation_series.R",local=TRUE)
source(file="./../../QML_Variance.R",local=TRUE)

theta_1a <- c(0.0001, 0.12, 0.85)
theta_2a <- c(0.0001, 0.32, 0.65)
theta_3a <- c(0.0001, 0.08, 0.82)
epsilon <- simulation_rendements(2000, theta_1a)[1500:2000]
epsilon_chgt <- simulation_rendements_avec_changement_GARCH(2000, theta_1a, theta_2a, cut=0.8)[1500:2000]

eps1 <- simulation_rendements(2000, theta_1a)
#eps2 <- simulation_rendements(2000, theta_2a)
#eps3 <- simulation_rendements(2000, theta_3a)
#plot(eps1,type='l')
#lines(0:1999, eps2, col='blue')
#lines(0:1999, eps3, col="red")

#sigma2_1 <- simu_sigma2(eps1**2, theta_1a)
#sigma2_2 <- simu_sigma2(eps2**2, theta_2a)
#sigma2_3 <- simu_sigma2(eps3**2, theta_3a)

#plot(sigma2_1, type="l")
#lines(sigma2_2,col="blue")
#lines(sigma2_3, col="red")

sigma2_1 <- simu_sigma2(eps1**2, theta_1a)
sigma2_1_2a <- simu_sigma2(eps1**2, theta_2a)
sigma2_1_3a <- simu_sigma2(eps1**2, theta_3a)
#plot(sigma2_1, type='l', col="red")
#lines(sigma2_1_2a, col="blue")
#lines(sigma2_1_3a, col="green")

#GGploting2
theta_1a <- c(0.0001, 0.12, 0.85)
theta_2a <- c(0.0001, 0.32, 0.65)
theta_3a <- c(0.0001, 0.08, 0.82)
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
plott <- plott + geom_line(aes(y=False_sigma2_a),linetype=2, color="red")
plott

plott2 <- ggplot(data=dtfrm, aes(x=Date))
plott2 <- plott2 + geom_line(aes(y=True_sigma2), linetype=2, color="blue")
plott2 <- plott2 + geom_line(aes(y=False_sigma2_b),linetype=2, color="red")
plott2
#Stop ggploting2

epsilon <- epsilon_chgt

simu_sigma2_eps <- simu_sigma2(epsilon**2, theta_1a)

mat_sigma <- as.matrix(simu_sigma2_eps, ncol=1)
indice <- as.matrix(1:length(mat_sigma), ncol=1)
mat_sigma <- cbind(indice, mat_sigma)

sigma_2 <- as.data.frame(mat_sigma)
colnames(sigma_2) <- c("Date","Sigma2")
plotting <- ggplot(data = sigma_2, aes(x = Date, y = Sigma2))
plotting <- plotting + geom_line(color="blue")

asymptote <- cbind(indice, as.matrix(rep(theta_1a[1]/(1-theta_1a[2]-theta_1a[3]), length(indice)), ncol=1))
asymptote <- as.data.frame(asymptote)
colnames(asymptote) <- c("Date","Asymptote")
plotting <- plotting + geom_line(data = asymptote, aes(x = Date, y = Asymptote))
plotting
