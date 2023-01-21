library(dplyr)

data = read.csv("C:/Users/maeld/OneDrive/Documents/Ensae/2A/Statapp/donnees/DAX.csv")
source(file = "C:/Users/maeld/OneDrive/Documents/GitHub/Garch_Model_Assets/data_preparation.R",local= TRUE)
data <- transform_csv(data)

source(file = "C:/Users/maeld/OneDrive/Documents/GitHub/Garch_Model_Assets/QML_Mael.R",local= TRUE)
plot(data$rendement2)
eps2 = data$rendement2[2:length(data$rendement2)]
eps2 = eps2/sd(eps2)
n = length(eps2)
theta = QML(eps2[1:200])$par #pb!
print(theta)


init = mean(eps2[1:200])
sigma2 = simu_sigma2(init,eps2[1:n],theta)

Z = c()
for(i in 201:n){
  if((eps2[i]/sigma2[i])<1.96**2) {Z[i-200]=1}
  else {Z[i-200]=0}  }

plot(Z)

#test d'adéquation du khi2 sur la binomiale
inside = length(Z[Z==1])
outside = length(Z[Z==0])
res = c(outside,inside)
proba = c(0.05,0.95)
chisq.test(res,p=proba)
  
