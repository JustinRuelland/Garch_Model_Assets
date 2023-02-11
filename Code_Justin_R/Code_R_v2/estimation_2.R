#----------------- Estimation horizon 2 --------------------
options( "digits"=7, "scipen"=10) 
rm(list=ls())
library(signal)
library(dplyr)
library(ggplot2)
library(lubridate)

source(file = "../../data_preparation.R",local= TRUE) 
source(file = "../../QML_Variance.R",local= TRUE)

#Import du cac40
data = read.csv("../../CAC40_15_19.csv")
data <- transform_csv(data)
plot_series_temp(data)


#Extraction epsilon square
epsilon_square <- data$rendement2
epsilon <- data$rendement

#Estimation des paramètres du modèle Garch fit sur 80% de la serie
n_eps <- length(epsilon)
threshold_train <- floor(0.8*n_eps)
epsilon_train <- epsilon[1:threshold_train]

theta_hat <- QML(epsilon_train**2)

epsilon_test <- epsilon[(threshold_train+1):n_eps]

# Serie a traiter (à mettre sous forme de fonction et boucle)
epsilon_traited <- c(epsilon_train, epsilon_test[1:1])
epsilon_square_traited <- epsilon_traited**2

# Obtention des sigmas_square et des eta_square
sigma_square_traited <- simu_sigma2(epsilon_traited**2, theta_hat)
sigma_traited <- sqrt(sigma_square_traited)
eta_traited <- epsilon_traited/sigma_traited

n_long <- threshold_train + 1

boostraped_estimation <- function(epsilon_traited, sigma_traited, eta_traited, n_long, n_path, n_day){
  simulated_garch <- matrix(data = 0, nrow = n_path, ncol =100)
  for(j in 1:n_path){
    #Bootstrap for n etas
    eta_bootstraped <- sample(eta_traited, size = n_day, replace = TRUE)
    
    #Following for n sigmas
    sigmas_square <- c(sigma_square_traited)
    etas <- c(eta_traited, eta_bootstraped)
    eps_simulated <- c(epsilon_traited)
    n_eps_traited <- length(sigma_traited)
    for(i in 1:n_day){
      k <- n_eps_traited + i
      sigmas_square[k] = theta_hat[1] + theta_hat[2]*(sigmas_square[k-1]*etas[k-1]**2) + theta_hat[3]*sigmas_square[k-1]
      eps_simulated[k] = sqrt(sigmas_square[k]) * etas[k]
    }
    cut <- length(eps_simulated) - 99
    simulated_garch[j,] = eps_simulated[cut:length(eps_simulated)]
  }
  return(simulated_garch)
}

path <- 300
cutting <- 10
cut <- n_long-77
cut_end <- n_long+22
matrix_sim <- cbind(t(boostraped_estimation(epsilon_traited, sigma_traited, eta_traited, n_long, path, 70)), cut:cut_end)
simulation <- as.data.frame(matrix_sim)
colnames(simulation)[path+1] <- "Date" 
t <- "Date"

p <- ggplot(data=simulation)
rainbow_p <- rainbow(path*100, s=.6, v=.9)[sample(1:(path*100), path)]
ic_garch <- rep(0.016, 30)
for(i in 1:path){
  col <- paste("V", as.character(i), sep="")
  p <- p + geom_line(aes(x=Date, y= .data[[col]]), color=rainbow_p[i])
  #p <- p + geom_line(aes(x=Date, y= .data[[col]]), color="black")
  
}
for(i in 1:70){
  ic_garch <- c(ic_garch, quantile(matrix_sim[i,1:300], probs=0.95))
  ic_garch_2 <- c(ic_garch, quantile(matrix_sim[i,1:300], probs=0.05))
}
#p <- ggplot(data=simulation) + geom_line(aes(x=.data[[t]], y=V1), color="blue") + geom_line(aes(x=Date, y=V2), color="blue") + geom_line(aes(x=Date, y=V2), color="blue")
#p <- p + geom_line(aes(x=Date, y= ))
#Bootstrap for 2 etas

p <- p + geom_line(aes(x=Date, y=0.016), color="black")
p <- p + geom_line(aes(x=Date, y=-0.018), color="black")
c2 <- cut+48
c2_End <- cut_end + 48
p <- p + geom_line(aes(x=Date, y=epsilon[c2:c2_End]), color="black")
#p <- p + geom_lin
p <- p + geom_line(aes(x=Date, y=ic_garch[1:100]), color="brown")
p <- p + geom_line(aes(x=Date, y=ic_garch_2[1:100]), color="brown")
p

eta_bootstraped <- sample(eta_traited, size = 2, replace = TRUE)

#Following for 2 sigmas
sigmas_square <- c(sigma_square_traited)
etas <- c(eta_traited, eta_bootstraped)
eps_simulated <- c(epsilon_traited)
n_eps_traited <- length(sigma_traited)
for(i in 1:2){
  k <- n_eps_traited + i
  sigmas_square[k] = theta_hat[1] + theta_hat[2]*(sigmas_square[k-1]*etas[k-1]**2) + theta_hat[3]*sigmas_square[k-1]
  eps_simulated[k] = sqrt(sigmas_square[k]) * etas[k]
}



