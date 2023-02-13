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

trajectories_bootstrapped <- function(sigma_init, eps_init, etas_estimated, n_path, n_day, theta){
  matrix_simulation <- matrix(data = 0, nrow = n_path, ncol = n_day)
  for(i in 1:n_path){
    etas <- bootstrap_etas(etas_estimated, n_day)
    
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

day_simulated <- 500
day_displaying <- 1000 # take into account the simulated days
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

#OOOF
matrix_simulation = matrix_trajectories
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
#p <- p + geom_line(aes(x=df[["Index"]], y=df[["V1"]]), color=rainbow_p[1])
#p <- p + geom_line(aes(x=df[["Index"]], y=df[["V2"]]), color=rainbow_p[2])
#p <- p + geom_line(aes(x=df[["Index"]], y=df[["V3"]]), color=rainbow_p[3])
#p <- p + geom_line(aes(x=df[["Index"]], y=df[["V4"]]), color=rainbow_p[4])
#p <- p + geom_line(aes(x=df[["Index"]], y=df[["V5"]]), color=rainbow_p[5])
#p <- p + geom_line(aes(x=df[["Index"]], y=df[["V6"]]), color=rainbow_p[6])
#p <- p + geom_line(aes(x=df[["Index"]], y=df[["V7"]]), color=rainbow_p[7])
#p <- p + geom_line(aes(x=df[["Index"]], y=df[["V8"]]), color=rainbow_p[8])
#p <- p + geom_line(aes(x=df[["Index"]], y=df[["V9"]]), color=rainbow_p[9])
#m <- p + geom_line(aes(x=df[["Index"]], y=df[["V1"]]), color=rainbow_p[1])
for(i in 1:n_path){
  col <- paste("V", as.character(i), sep = "")
  p <- p + geom_line(aes_(x=df[["Index"]], y=df[[col]]), color=rainbow_p[i])
}
p <- p + geom_line(aes_(x = df_true[["Date"]], y = df_true[["V1"]]), color="black")
p
#OOOF

estimation_day <- threshold_train+1
end_day <- estimation_day + day_simulated - 1
start_day <- estimation_day - (day_displaying - day_simulated)
epsilon_initial <- epsilon[start_day:end_day]

displaying_trajectories(epsilon_initial = epsilon_initial, matrix_simulation = matrix_trajectories, day_displaying = day_displaying)



# Studying of the trajectories

matrix_trajectories <- matrix(data=0, nrow=path_simulated, ncol = day_simulated)
matrix_trajectories <- trajectories_bootstrapped(sigma_init = sigma_init, eps_init = eps_init, etas_estimated = eta_traited, n_path = path_simulated, n_day = day_simulated, theta = theta_hat)

#---Moyenne---
Mean_along_sim <- matrix(data=NA, nrow = 1, ncol = day_simulated)

#---IC 95%---
IC_along_sim_up <- matrix(data=NA, nrow = 1, ncol = day_simulated)
IC_along_sim_down <- matrix(data=NA, nrow = 1, ncol = day_simulated)

for(i in 1:day_simulated){
  Mean_along_sim[1,i] <- mean(matrix_trajectories[,i])
  IC_along_sim_up[1,i] <- quantile(matrix_trajectories[,i], probs = 0.975)
  IC_along_sim_down[1,i] <- quantile(matrix_trajectories[,i], probs = 0.025)
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
p2 <- p2 + geom_line(aes(x=df_ic_up_naif[["Date"]], y=df_ic_up_naif[["V1"]]), color="orange")
p2 <- p2 + geom_line(aes(x=df_ic_down_naif[["Date"]], y=df_ic_down_naif[["V1"]]), color="orange")
p2 <- p2 + ggtitle("Estimation APPL horizon 500 jours par modèle GARCH") + xlab("Date") + ylab("Rendement")
p2





# Old version of the code


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



