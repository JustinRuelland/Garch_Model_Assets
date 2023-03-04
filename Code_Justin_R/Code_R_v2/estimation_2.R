#----------------- Estimation horizon 2 --------------------
options( "digits"=15, "scipen"=10) 
rm(list=ls())
library(signal)
library(dplyr)
library(ggplot2)
library(lubridate)
setwd("/Users/justinr/Documents/2A/Stat App/Github/Garch_Model_Assets/Code_Justin_R/Code_R_v2")
source(file = "../../data_preparation.R",local= TRUE) 
source(file = "../../QML_Variance.R",local= TRUE)

#Import du cac40
data = read.csv("../Yahoo_import/AAPL.csv")
data <- transform_csv(data)
plot_series_temp(data)


#Extraction epsilon square
epsilon_square <- data$rendement2
epsilon <- data$rendement

#Estimation des paramètres du modèle Garch fit sur 80% de la serie
n_eps <- length(epsilon)
threshold_train <- floor(0.5*n_eps)
epsilon_train <- epsilon[1:threshold_train]

theta_hat <- QML(epsilon_train**2)

epsilon_test <- epsilon[(threshold_train+1):n_eps]

# Serie a traiter (à mettre sous forme de fonction et boucle)
epsilon_traited <- c(epsilon_train, epsilon_test[1])
epsilon_square_traited <- epsilon_traited**2

# Obtention des sigmas_square et des eta_square
sigma_square_traited <- simu_sigma2(epsilon_traited**2, theta_hat)
sigma_traited <- sqrt(sigma_square_traited)
eta_traited <- epsilon_traited/sigma_traited

# Fonctions d'obtention des etas
bootstrap_etas <- function(etas, n){
  return(sample(x = etas, size = n, replace = TRUE))
}

# Fitting de densité
hist(eta_traited, freq = F)
lines(density(eta_traited), col = "red")
xxx <- seq(from = min(hist(eta_traited, plot = F)$breaks), to=max(hist(eta_traited, plot = F)$breaks), by=0.1)
lines(x = xxx, y = dnorm(x = xxx, mean = mean(eta_traited), sd = sd(eta_traited)), col = "blue")



fitting_etas <- function(etas, n){
  u <- runif(n)
  return(quantile(etas, u))
}



# Simulation des bras inconnus
trajectory <- function(sigma_init, eps_init, etas, theta){
  n <- length(etas)
  eps <- 1:n

  sigma_square <- theta[1] + theta[2]*(eps_init**2) + theta[3]*(sigma_init**2)
  
  for(i in 1:n){
    eps[i] <- sqrt(sigma_square)*etas[i]
    
    # Computation of the sigma_square_{i+1}
    sigma_square <- theta[1] + theta[2]*(eps[i]**2) + theta[3]*(sigma_square)
  }
  
  return(eps)
}

trajectories_bootstrapped <- function(sigma_init, eps_init, etas_estimated, n_path, n_day, theta, bootstrap=TRUE){
  matrix_simulation <- matrix(data = 0, nrow = n_path, ncol = n_day)
  for(i in 1:n_path){
    if(bootstrap){
      etas <- bootstrap_etas(etas_estimated, n_day)
    }else{
      etas <- fitting_etas(etas_estimated, n_day)
    }
    
    
    eps <- 1:n_day
    sigma_square <- theta[1] + theta[2]*(eps_init**2) + theta[3]*(sigma_init**2)
    
    for(j in 1:n_day){
      eps[j] <- sqrt(sigma_square)*etas[j]
  
      sigma_square <- theta[1] + theta[2]*(eps[j]**2) + theta[3]*(sigma_square)
    }
    
    matrix_simulation[i,] <- eps
  }
  return(matrix_simulation)
}

# Display of the trajectories

day_simulated <- 30
day_displaying <- 50 # take into account the simulated days
path_simulated <- 3000

sigma_init <- sigma_traited[length(sigma_traited)]
eps_init <- epsilon_traited[length(epsilon_traited)]

matrix_trajectories <- matrix(data=0, nrow=path_simulated, ncol = day_simulated)
matrix_trajectories <- trajectories_bootstrapped(sigma_init = sigma_init, eps_init = eps_init, etas_estimated = eta_traited, n_path = path_simulated, n_day = day_simulated, theta = theta_hat)

displaying_trajectories <- function(epsilon_initial, matrix_simulation, day_displaying){
  start <- day_displaying - length(matrix_simulation[1,]) + 1
  matrix_displaying <- matrix(data = NA, nrow=path_simulated, ncol = (day_displaying - day_simulated))
  matrix_displaying[,(start-1)] <- rep(epsilon_initial[(start-1)], length(matrix_simulation[,1]))
  matrix_displaying <- cbind(matrix_displaying, matrix_simulation)
  matrix_displaying <- cbind(t(matrix_displaying), 1:day_displaying)
  
  df <- as.data.frame(matrix_displaying)
  n_path <- length(matrix_simulation[,1])
  colnames(df)[n_path+1] <- "Index"
  
  df_true <- as.data.frame(cbind(as.matrix(epsilon_initial, ncol=1), 1:day_displaying))
  colnames(df_true)[2] <- "Date"
  
  rainbow_p <- rainbow(n_path, s=.6, v=.9)[sample(1:n_path, n_path)]
  p <- ggplot(data = df_true)
  p <- p + geom_line(aes(x = Date, y = V1), color="black")
  for(i in 1:n_path){
    col <- paste("V", as.character(i), sep = "")
    p <- p + geom_line(aes_(x=df[["Index"]], y=df[[col]]), color=rainbow_p[i])
  }
  p <- p + geom_line(aes_(x = df_true[["Date"]], y = df_true[["V1"]]), color="black")
  p
}


estimation_day <- threshold_train+1
end_day <- estimation_day + day_simulated - 1
start_day <- estimation_day - (day_displaying - day_simulated)
epsilon_initial <- epsilon[start_day:end_day]

displaying_trajectories(epsilon_initial = epsilon_initial, matrix_simulation = matrix_trajectories, day_displaying = day_displaying)



# Studying of the trajectories

matrix_trajectories <- matrix(data=0, nrow=path_simulated, ncol = day_simulated)
matrix_trajectories <- trajectories_bootstrapped(sigma_init = sigma_init, eps_init = eps_init, etas_estimated = eta_traited, n_path = path_simulated, n_day = day_simulated, theta = theta_hat)
matrix_trajectories_notbootstrapped <- trajectories_bootstrapped(sigma_init = sigma_init, eps_init = eps_init, etas_estimated = eta_traited, n_path = path_simulated, n_day = day_simulated, theta = theta_hat, bootstrap = FALSE)

#---Moyenne---
Mean_along_sim <- matrix(data=NA, nrow = 1, ncol = day_simulated)
Mean_along_sim_nb <- matrix(data=NA, nrow = 1, ncol = day_simulated)

#---IC 95%---
IC_along_sim_up <- matrix(data=NA, nrow = 1, ncol = day_simulated)
IC_along_sim_down <- matrix(data=NA, nrow = 1, ncol = day_simulated)

for(i in 1:day_simulated){
  Mean_along_sim[1,i] <- mean(matrix_trajectories[,i])
  IC_along_sim_up[1,i] <- quantile(matrix_trajectories[,i], probs = 0.975)
  IC_along_sim_down[1,i] <- quantile(matrix_trajectories[,i], probs = 0.025)
}

#---IC 95%--- not boostrapped
IC_along_sim_up_nb <- matrix(data=NA, nrow = 1, ncol = day_simulated)
IC_along_sim_down_nb <- matrix(data=NA, nrow = 1, ncol = day_simulated)

for(i in 1:day_simulated){
  Mean_along_sim_nb[1,i] <- mean(matrix_trajectories_notbootstrapped[,i])
  IC_along_sim_up_nb[1,i] <- quantile(matrix_trajectories_notbootstrapped[,i], probs = 0.975)
  IC_along_sim_down_nb[1,i] <- quantile(matrix_trajectories_notbootstrapped[,i], probs = 0.025)
}

#---IC 95% naif---
IC_naif_up <- matrix(data = NA, nrow=1, ncol = day_simulated)
quantile_up <- quantile(epsilon_initial[1:(day_displaying - day_simulated)], probs = 0.975)
IC_naif_up[1,] <- quantile_up

IC_naif_down <- matrix(data = NA, nrow=1, ncol = day_simulated)
quantile_down <- quantile(epsilon_initial[1:(day_displaying - day_simulated)], probs = 0.025)
IC_naif_down[1,] <- quantile_down

# Displaying these indicators
df_true <- as.data.frame(cbind(as.matrix(epsilon_initial, ncol=1), 1:day_displaying))
colnames(df_true)[2] <- "Date"

df_mean <- matrix(data = NA, nrow=1, ncol = (day_displaying - day_simulated))
df_mean<- cbind(df_mean, Mean_along_sim)
df_mean <- cbind(t(df_mean), 1:day_displaying)
df_mean <- as.data.frame(df_mean)
colnames(df_mean)[2] <- "Date"

df_ic_up <- matrix(data = NA, nrow=1, ncol = (day_displaying - day_simulated))
df_ic_up<- cbind(df_ic_up, IC_along_sim_up)
df_ic_up <- cbind(t(df_ic_up), 1:day_displaying)
df_ic_up <- as.data.frame(df_ic_up)
colnames(df_ic_up)[2] <- "Date"

df_ic_down <- matrix(data = NA, nrow=1, ncol = (day_displaying - day_simulated))
df_ic_down<- cbind(df_ic_down, IC_along_sim_down)
df_ic_down <- cbind(t(df_ic_down), 1:day_displaying)
df_ic_down <- as.data.frame(df_ic_down)
colnames(df_ic_down)[2] <- "Date"

#not bootstrappe
df_mean_nb <- matrix(data = NA, nrow=1, ncol = (day_displaying - day_simulated))
df_mean_nb<- cbind(df_mean_nb, Mean_along_sim_nb)
df_mean_nb <- cbind(t(df_mean_nb), 1:day_displaying)
df_mean_nb <- as.data.frame(df_mean_nb)
colnames(df_mean_nb)[2] <- "Date"

df_ic_up_nb <- matrix(data = NA, nrow=1, ncol = (day_displaying - day_simulated))
df_ic_up_nb<- cbind(df_ic_up_nb, IC_along_sim_up_nb)
df_ic_up_nb <- cbind(t(df_ic_up_nb), 1:day_displaying)
df_ic_up_nb <- as.data.frame(df_ic_up_nb)
colnames(df_ic_up_nb)[2] <- "Date"

df_ic_down_nb <- matrix(data = NA, nrow=1, ncol = (day_displaying - day_simulated))
df_ic_down_nb<- cbind(df_ic_down_nb, IC_along_sim_down_nb)
df_ic_down_nb <- cbind(t(df_ic_down_nb), 1:day_displaying)
df_ic_down_nb <- as.data.frame(df_ic_down_nb)
colnames(df_ic_down_nb)[2] <- "Date"
#---

df_ic_up_naif <- matrix(data = NA, nrow=1, ncol = (day_displaying - day_simulated))
df_ic_up_naif<- cbind(df_ic_up_naif, IC_naif_up)
df_ic_up_naif <- cbind(t(df_ic_up_naif), 1:day_displaying)
df_ic_up_naif <- as.data.frame(df_ic_up_naif)
colnames(df_ic_up_naif)[2] <- "Date"

df_ic_down_naif <- matrix(data = NA, nrow=1, ncol = (day_displaying - day_simulated))
df_ic_down_naif<- cbind(df_ic_down_naif, IC_naif_down)
df_ic_down_naif <- cbind(t(df_ic_down_naif), 1:day_displaying)
df_ic_down_naif <- as.data.frame(df_ic_down_naif)
colnames(df_ic_down_naif)[2] <- "Date"

p2 <- ggplot(data = df_true) + geom_line(aes(x = Date, y = V1), color="black")
p2 <- p2 + geom_line(aes(x=df_mean[["Date"]], y=df_mean[["V1"]]), color="blue")
p2 <- p2 + geom_line(aes(x=df_ic_up[["Date"]], y=df_ic_up[["V1"]]), color="red")
p2 <- p2 + geom_line(aes(x=df_ic_down[["Date"]], y=df_ic_down[["V1"]]), color="red")
#NB
p2 <- p2 + geom_line(aes(x=df_mean_nb[["Date"]], y=df_mean_nb[["V1"]]), color="green")
p2 <- p2 + geom_line(aes(x=df_ic_up_nb[["Date"]], y=df_ic_up_nb[["V1"]]), color="purple")
p2 <- p2 + geom_line(aes(x=df_ic_down_nb[["Date"]], y=df_ic_down_nb[["V1"]]), color="purple")
#--
p2 <- p2 + geom_line(aes(x=df_ic_up_naif[["Date"]], y=df_ic_up_naif[["V1"]]), color="orange")
p2 <- p2 + geom_line(aes(x=df_ic_down_naif[["Date"]], y=df_ic_down_naif[["V1"]]), color="orange")
p2 <- p2 + ggtitle("Estimation APPL horizon 500 jours par modèle GARCH") + xlab("Date") + ylab("Rendement")
p2


# Visualisation simulation par méthode d'inversion
density_fitted <- density(eta_traited)
plot(density_fitted, col="red")
eta_sim <- runif(629)
for(i in 1:629){
  eta_sim[i] <- quantile(eta_traited, eta_sim[i])
}
lines(density(eta_sim), col="blue")
