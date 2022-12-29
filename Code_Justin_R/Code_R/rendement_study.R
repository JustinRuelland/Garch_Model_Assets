library(ggplot2)
library(dplyr)
data <- read.csv2("./../Yahoo_import/FTT-USD.csv")
data_2 <- read.csv2("./../Yahoo_import/AAPL.csv")

data_2 <- read.csv2("./../Yahoo_import/GOOG.csv")

data <- data$Date.Open.High.Low.Close.Adj.Close.Volume
data_2 <- data_2$Date.Open.High.Low.Close.Adj.Close.Volume

x_quotes <- c()
x_quotes_2 <- c()

for(i in data){
  x_quotes <- c(x_quotes, as.numeric(strsplit(i, ",")[[1]][2]))
}

for(i in data_2){
  x_quotes_2 <- c(x_quotes_2, as.numeric(strsplit(i, ",")[[1]][2]))
}

plot(x_quotes, type = "l")
plot(x_quotes_2, type = "l")

## Etude de AAPL sur les 500 derniers jours -----

x_quotes_2_cur <- x_quotes_2[500:1258]
plot(x_quotes_2_cur, type='l', main = "Prix de l'actif DAXX50 sur les 500 derniers jours")

## Etude des rendements
epsilon_2 <- log(x_quotes[2:200]/x_quotes[1:199])
epsilon = log(x_quotes_2_cur[2:500]/x_quotes_2_cur[1:499])

plot(epsilon, type="l", col="blue", main = "Rendements en log de l'actif AAPL")
plot(epsilon_2, type="l")

epsilon <- epsilon**2
## Etude des moments des rendements
epsilon <- serie
covar_biased <- function(vec1, vec2){
  n <- length(vec1)
  return((1/n)*sum((vec1 - mean(vec1))*(vec2 - mean(vec2))))
}

mu <- mean(epsilon)
Variance <- covar_biased(epsilon, epsilon)

N <- length(epsilon)

cov_hat <- function(h){
  rate <- N-h
  dem_h <- h+1
  eps_centered <- epsilon 
  eps_t <-  eps_centered[1:rate]
  eps_th <- eps_centered[dem_h:N]
  return(covar_biased(eps_t, eps_th))
}

corr_hat <- function(h){
  return(cov_hat(h)/Variance)
}


corr_hat <- Vectorize(FUN = corr_hat, vectorize.args = "h")
cov_hat <- Vectorize(FUN = cov_hat, "h")

corr_hat_2 <- corr_hat(0:498)

plot(1:30, corr_hat_2[1:30])
lines(x = 0:100, y=rep(0,101), )
#lines(x = 0:100, y = rep(1.96/sqrt(498), 101), col="blue")
#lines(x = 0:100, y = rep(-1.96/sqrt(498), 101), col="blue")


x <- rnorm(100000)
x <- x**2

ligne <- c()
for(a in 1:400){
  alpha <- (a-1)/100
  for(b in 1:200){
    beta <- (b-1)/100
    y <- mean(log(alpha*x+beta))
    y_2 <- var(log(alpha*x+beta))
    if(y >= 0){
      ligne <- c(ligne, beta)
      break
    }
  }
}

xx <- seq(from=0, to=3.99, by=0.01)
ggplot(as.data.frame(ligne)) + geom_line(aes(xx, y=ligne)) + ylim(0,1.5) + xlim(0,4) + xlab("alpha") + ylab("beta") + ggtitle("Condition de stricte stationnarité") + theme(plot.title=element_text(hjust=0.5))



#---------------------- Echantillon_norm epsilon ---------------------------
rm(list=ls())

omega_0 <- 0.5
alpha_0 <- 0.12
beta_0 <- 0.85

taille_relative_echantillon <- 3000
Nb_echantillon <- 5000

a_lin <- function(x, alpha, beta){
  return(alpha * x + beta)
}

a_lin <- Vectorize(FUN = a_lin, vectorize.args = "x")

echantillon_eta_square_1000 <- matrix(data = rnorm(n= taille_relative_echantillon*Nb_echantillon)**2, nrow = Nb_echantillon)

echant_eta_square_1000 <- as.list(as.data.frame(t(echantillon_eta_square_1000)))

#---------------------- Computation of sigma square ------------------------
a = 0

splitting <- 1000

sigma_square_0 <- function(eta, omega, alpha, beta){
  vec <- a_lin(eta[1:splitting], alpha = alpha, beta = beta)
  somme <- 0
  for(i in 1:splitting){
    indice <- splitting - i + 1
    indice_end <- splitting - 1
    somme = somme + prod(vec[indice:indice_end])
  }
  a <<- a+1
  print(a)
  return((1 + somme)*omega)
}

sigma_square_0 <- Vectorize(FUN = sigma_square_0, vectorize.args = "eta")

sigma_sq_0 <- sigma_square_0(echant_eta_square_1000, omega_0, alpha_0, beta_0)

#---------------------- Computation of epsilon square ----------------------
compteur_2 <- 0

serie_rendement <- function(eta, sigma_init, omega, alpha, beta){
  epsilon <- eta[splitting]*sigma_init
  last_eps <- epsilon
  last_sigm <- sigma_init
  starting <- splitting + 1
  for(i in starting:taille_relative_echantillon){
    sigma <- omega + alpha*last_eps + beta*last_sigm
    last_sigm <- sigma
    eps_curr <- sigma*eta[i]
    last_eps <- eps_curr
    epsilon <- c(epsilon, last_eps)
  }
  compteur_2 <<- compteur_2 + 1
  print(compteur_2)
  return(epsilon)
}

#serie_rendement <- Vectorize(FUN = serie_rendement, vectorize.args = c("eta", "sigma_init"))

serie <- serie_rendement(eta = echant_eta_square_1000[[5]], sigma_init = sigma_sq_0[5], omega = omega_0, alpha = alpha_0, beta = beta_0)


plot(1:2001, serie, type='l')


#---------------------- Methode du quasi-maximum de vraisemblance ----------


creat_vec <- function(scalaire, niveau, taille){
  size <- taille - niveau
  return( c(rep(1,niveau), rep(scalaire, size)))
} #fonction outil de la fonction expansion log

expansion_log <- function(suites){
  taille <- length(suites)
  lim <- taille - 1
  scalaire <- suites[1]
  for(i in 1:lim){
    suites <- suites * creat_vec(scalaire, i, taille);
  }
  return(c(1,suites))
} #fonction de calcul du terme beta mis en puissance


estimation_sigma_initial <- function(past, omega, alpha, beta){
  constante <- omega/(1-beta)
  eviction <- past * alpha
  taille <- length(past)
  eviction = eviction * expansion_log(rep(beta, taille-1))
  return(sum(eviction)+constante)
}

serie_sigma <- function(sigma_init, serie, omega, alpha, beta){
  sigma <- sigma_init
  taille <- length(serie) - 1
  last_sigm <- sigma
  for(i in 1:taille){
    sigma_curr <- omega + alpha*serie[i] + beta*last_sigm
    last_sigm <- sigma_curr
    sigma <- c(sigma, sigma_curr)
  }
  return(sigma)
}

quasi_max_likelihood <- function(omega, alpha, beta, serie){
  taille <- length(serie);
  echantillon_cible <- floor((1/6)*taille);
  
  #Calcul de la serie des sigma a partir des epsilons
  sigma_init <- estimation_sigma_initial(past = serie[1:echantillon_cible], omega = omega, alpha = alpha, beta = beta)
  indice_1 <- echantillon_cible + 1
  serie_study <- serie[indice_1:taille]
  sigma <- serie_sigma(sigma_init = sigma_init, serie = serie_study,  omega = omega, alpha = alpha, beta = beta)
  
  #Calcul des termes generaux de la serie
  vec <- log(sigma) + (serie_study)/sigma
  
  return(sum(vec))
}

quasi_max <- function(x, serie){
  return(quasi_max_likelihood(omega =  x[1], alpha =  x[2], beta =  x[3], serie = serie))
}

ui <- cbind(c(1,-1,0,0,0,0),c(0,0,1,0,-1,0),c(0,0,0,1,-1,-1))

ci <- c(0.01, -1, 0, 0, -1, -0.999)

test1 <- serie

optimisation <- constrOptim(theta = c(0.2,0.8,0.1), f = quasi_max, grad = NULL, ui = ui, ci = ci, serie = test1)

#---------------------- Analyse des paramètres estimés ---------------------

omega_hat <- c()
alpha_hat <- c()
beta_hat <- c()

for(i in 1:100){
  serie_test <- serie_rendement(eta = echant_eta_square_1000[[i]], sigma_init = sigma_sq_0[i], omega = omega_0, alpha = alpha_0, beta = beta_0)
  optimisation <- constrOptim(theta = c(0.5,0.4,0.4), f = quasi_max, grad = NULL, ui = ui, ci = ci, serie = serie_test)
  omega_hat <- c(omega_hat, optimisation$par[1])
  alpha_hat <- c(alpha_hat, optimisation$par[2])
  beta_hat <- c(beta_hat, optimisation$par[3])
}

omega <- omega_hat - omega_0
alpha <- alpha_hat - alpha_0
beta <- beta_hat - beta_0

matrice <- as.matrix(cbind(omega, alpha, beta))
library(ggplot2)

p <- ggplot(as.data.frame(matrice)) + geom_boxplot()
