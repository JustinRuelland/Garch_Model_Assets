rm(list=ls())
setwd(dir = "/Users/justinr/Documents/2A/Stat App/Github/Garch_Model_Assets/Code_Justin_R/Code_R_v2")
library(ggplot2)
library(dplyr)
library(forecast)
library(keras)
library(tensorflow)

source(file="./../../simulation_series.R",local=TRUE)
source(file="./../../QML_Variance.R",local=TRUE)
source(file="./../../data_preparation.R",local=TRUE)

# Regression linéaire 

# Import series financieres reelles
data_real <- read.csv("./../../^GDAXI.csv") #fichier csv de Yahoo finance
data_real <- transform_csv(data_real)

eps_2 <- simulation_rendements(n = 2000, theta = c(0.00001,0.1,0.85))**2

QML(eps_2)

eps_square <- eps_2
eps_square <- data_real$rendement2

# Dropping na
eps_square <- eps_square[!is.na(eps_square)]

# Visualisation des correlograms pour determiner les variables d'interets
autocorrelations_korrigieren <-function(data){
  #acf(data$rendement,type='correlation', na.action=na.pass, lag.max = 0,plot=TRUE)
  acf(data$rendement2,type='correlation', na.action=na.pass, lag.max = 100,plot=TRUE)
}


autocorrelations_korrigieren(data_real)
# On identifie une baisse de correlation au dela de 40

#------- Creation du database -------  

# Tools functions for naming columns

colname_erstellung <- function(lag.max){
  naming <- as.character(1:lag.max)
  for(i in 1:lag.max){
    naming[i] <- paste(naming[i], "lag")
  }
  return(naming)
}

colname_erstellung_prediction <- function(lag.max){
  naming <- "Intercept"
  return(c(naming, colname_erstellung(lag.max)))
}

# Tools functions for creating df

database_erstellung <- function(serie, lag.max){
  taille <- length(serie)
  nb_ligne <- taille - lag.max
  data_matrix_X <- matrix(NA, ncol = lag.max, nrow = nb_ligne)
  data_matrix_Y <- matrix(NA, ncol = 1, nrow = nb_ligne)
  debut <- lag.max + 1
  for(i in 1:nb_ligne){
    current <- i + lag.max
    data_matrix_Y[i] <- serie[current]
    first <- current - 1
    last <- current - lag.max
    data_matrix_X[i,] <- serie[first:last]
  }
  colnames(data_matrix_Y) <- "Y"
  colnames(data_matrix_X) <- colname_erstellung(lag.max)
  
  df <- as.data.frame(cbind(data_matrix_Y, data_matrix_X))
  return(df)
}

database_erstellung_prediction <- function(serie, lag.max){
  taille <- length(serie)
  nb_ligne <- taille - lag.max + 1
  nb_col <- lag.max + 1
  data_matrix_X <- matrix(NA, ncol = nb_col, nrow = nb_ligne)
  data_matrix_X[,1] <- rep(1, nb_ligne)
  for(i in 1:nb_ligne){
    current <- i + lag.max - 1
    first <- current
    last <- current - lag.max + 1
    data_matrix_X[i,2:nb_col] <- serie[first:last]
  }
  colnames(data_matrix_X) <- colname_erstellung_prediction(lag.max)
  
  df <- as.data.frame(data_matrix_X)
  return(df)
}

# Enregistrement design matrix
df_design <- database_erstellung(eps_square, 40)

# OLS testing

Y <- as.matrix(df_design$Y)
X <- as.matrix(df_design[colnames(df_design) != "Y"])
OLS_spurious <- lm(Y ~ X)

summary(OLS_spurious) 


#------- OLS prediction ------- 

OLS_prediction <- function(train_set, evaluate_set, lag.max){
  #DF Train Construction + Training
  df_design <- database_erstellung(train_set, lag.max)
  Y <- as.matrix(df_design$Y)
  X <- as.matrix(df_design[colnames(df_design) != "Y"])
  OLS_spurious <- lm(Y ~ X)
  
  #DF Test Construction 
  df_design_train <- database_erstellung_prediction(evaluate_set, lag.max)
  X_test <- as.matrix(df_design_train)
  
  #Prediction
  prediction <- X_test %*% OLS_spurious$coefficients
  
  return(prediction)
}


#------ ?Unit? Test OLS_prediction -----

# Test on the DAXX of the prediction + graph 
# /!\ The prediction at this state is only for
# horizon 1 and for the naive OLS /!\


df <- as.data.frame(cbind(1:length(eps_square),eps_square))
colnames(df) <- c("Date","ReturnSquare")

#cut for the test/train 
#the cut is use to train, but the test set contains some train parts that are used as variables/predictors

cut <- 0.75 #in proportion of the set length
lag.max <- 40

end_train <- floor(0.75*length(df$ReturnSquare))
start_test <- end_train + 1 - lag.max
end_test <- length(df$ReturnSquare)
real_start <- start_test + lag.max

eps_train <- df$ReturnSquare[1:end_train]
eps_test <- df$ReturnSquare[start_test:end_test]
eps_real <- df$ReturnSquare[real_start:end_test]

eps_predicted <- OLS_prediction(eps_train, eps_test, lag.max)

#We drop the last one predicted because we don't have the true one to compare

eps_predicted <- eps_predicted[1:(length(eps_predicted)-1)]

#Displaying comparison 

#Mean of the set to have a compare to a very naive
eps_mean <- rep(mean(eps_train), length(eps_real))

df_display <- as.data.frame(cbind(1:length(eps_real), eps_real, eps_predicted, eps_mean))
colnames(df_display) <- c("Dates", "Real", "OLS_Prediction", "Mean")

disp <- ggplot(data=df_display, aes(x = Dates)) + geom_line(aes(y=Real), color="blue") + geom_line(aes(y=OLS_Prediction), color='cyan') + geom_line(aes(y=Mean), color="red")  
disp

#Test DM vs Garch

source(file="./../../comparison_model.R",local=TRUE)

# GARCH prediction horizon 1

pred_h1_garch_modified <- function(eps2,cut){
  
  n = length(eps2)
  n_cut = floor(n*cut)
  
  theta =  QML(eps2[1:n_cut])
  pred = double(n-n_cut+1)
  init = theta[1]/(1-theta[2]-theta[3])
  #pred[1] = func_sigma2(n_cut,init,eps2[1:n_cut],theta)
  pred[1] <- simu_sigma2(eps2[1:n_cut],theta)[n_cut]
  
  l = n-n_cut+1
  for(i in 2:l){
    pred[i] = theta[1]+theta[2]*eps2[n_cut+i-2]+theta[3]*pred[i-1]
    if(is.na(pred[i])){
      cat("NA appear at",i)
    }
    }
  return(pred[2:l])
}


eps_garch <- pred_h1_garch_modified(eps_square, cut)


test_mariano(eps_garch, eps_predicted, eps_real, hor = 1)


RMSE <- function(estimation, truth){
  n <- length(estimation)
  if(n != length(truth)){
    return(-1)
  }
  return(sqrt(mean((estimation-truth)**2)))
}

# Display
df_display <- as.data.frame(cbind(1:length(eps_real), eps_real, eps_predicted, eps_mean, eps_garch))
colnames(df_display) <- c("Dates", "Real", "OLS_Prediction", "Mean", "GARCH")

disp <- ggplot(data=df_display, aes(x = Dates)) + geom_line(aes(y=Real), color="blue") + geom_line(aes(y=OLS_Prediction), color='cyan') + geom_line(aes(y=Mean), color="red") + geom_line(aes(y=GARCH), color="green")  
disp
# End of Display

RMSE_garch <- RMSE(eps_garch, eps_real)
RMSE_ols <- RMSE(eps_predicted, eps_real)

RMSE_garch
RMSE_ols


#OLS lasso pour la suite

# Deep Learning 

# Construction of a model of deep learning with Keras API, sequential NN => first layer with 20 units and relu activation => second layer 10 units and relu activation => last layer 1 units with relu

df_train <- database_erstellung(eps_train, 40)
X_train <- as.matrix(df_train[colnames(df_train) != "Y"])
Y_train <- as.matrix(df_train$Y)

df_test <- database_erstellung(eps_test, 40)
X_test <- as.matrix(df_test[colnames(df_test) != "Y"])
Y_test <- as.matrix(df_test$Y)

nn_model <- keras_model_sequential()
nn_model %>% 
  layer_dense(units = 256, activation = "relu", input_shape = c(40)) %>%
  layer_dropout(rate=0.8) %>%
  layer_dense(units = 128, activation = "relu") %>%
  layer_dropout(rate=0.8) %>%
  layer_dense(units = 20, activation = "relu") %>%
  layer_dropout(rate=0.5) %>%
  layer_dense(units=1, activation = "linear")

#nn_model %>%
#  layer_dense(units = 1, activation="linear", input_shape = c(40))
  
nn_model %>% compile(loss='mean_squared_error', optimizer = optimizer_adam(), metrics=c('mean_squared_error', 'accuracy'))

summary(nn_model)

nn_model %>% fit(X_train, Y_train, epochs=100, batch_size=128, validation_split=0.2)

#nn_model %>% evaluate(X_test, Y_test)

#nn_model %>% predict(as.matrix(t(X[1,])))

nn_prediction <- predict(object = nn_model, X_test)

test_mariano(nn_prediction, eps_predicted, eps_real, hor = 1)
test_mariano(nn_prediction, eps_garch, eps_real, hor = 1)

RMSE(nn_prediction, eps_real)
RMSE_ols
RMSE_garch
