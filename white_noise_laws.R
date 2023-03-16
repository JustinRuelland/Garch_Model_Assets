# Loi des eta, pour qu'on puisse mettre autre chose que "rnorm"

runif_normalisee <-function(n){ # Loi uniforme normalisée
  return(runif(n,-sqrt(3),sqrt(3)))
} 


rt_8_normalisee <-function(n){ # Loi de Student à 8 degrés de liberté normalisée
  return(rt(n,8)/sqrt(8/(8-2)))
} 

normalised_student<-function(n){
  df = 5 #degree of freedom
  res = rt(n,df)/(sqrt(df/(df-2)))
}
