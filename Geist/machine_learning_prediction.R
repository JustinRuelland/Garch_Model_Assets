source(file = "../data_preparation.R",local= TRUE) #syntaxe pour aller dans le réportoire mère du fichier ML_prediction.R

#install.packages("keras")
library(keras)

#install.packages("mhsmm")
library(mvtnorm)
library(mhsmm)


### Data import
data = read.csv("../^GDAXI.csv")
data = transform_csv(data)



p = ggplot(data) + geom_line(aes(x = data$Date,y=data$Prix))
p
