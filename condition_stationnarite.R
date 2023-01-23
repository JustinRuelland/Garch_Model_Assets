library(ggplot2)
library(signal)
#---------------------- Condition de stationnarité -------------------------------
### But : Représentations des (alpha,beta) vérifiant la condition de stationnarité

### Utilisation de la fonction :
# Mettre en paramètres la fonction de simulation et ses paramètres
# Par exemple "loi_eta = rnorm", et les paramètres 0 et 1

condition_stationnarite <- function(loi_eta){
  N = 400 #nombre de alphas
  M = 10**4 #nombre de betas (ne pas dépasser 10**4 a priori, pour temps de calculs raisonnable)
  nsimu= 10**5 #nombre de simulations
  
  a = seq(0,3,length.out=N)
  b = seq(0,1,length.out=M)
  y = array(0,dim=c(M,N)) #résultat final
  
  for(i in range(nsimu)){
    #on évalue tous les alphas et betas sur une simulation
    xx = array(loi_eta(N*M), dim=c(N,M)) 
    y = y + log(t(a*xx**2)+b)/nsimu 
  }
  
  #0 si espérance négative, 1 si positive
  y[y>0]=0
  y[y<0]=1
  y = apply(y,2,FUN=mean) # sur les betas
  
  # Filtration de la série (pour avoir une courbe lisse)
  y = sgolayfilt(y, p =2, n=25) #p = degré des polynômes, n = nombre de points pour chaque polynôme
   
  #Passage en dataframe (pour utilisation ggplot2)
  df = cbind(a,y)
  df = as.data.frame(df)
  
  return(df)
  
}

superposition_3graphiques<-function(df1,df2,df3){
  Graph_condi_satio = ggplot(data = df1)+
    geom_line(aes(x = a,y = y))+
    xlab("Alpha") + 
    ylab("Beta") + 
    geom_line(data = df2,aes(x=a,y=y))+
    geom_line(data = df3,aes(x=a,y=y))
  
  plot(Graph_condi_satio)
  
}

# Test de la fonction
#condition_stationnarite(rnorm)
runif_normalisee <-function(n){return(runif(n,-sqrt(3),sqrt(3)))}
rt_8_normalisee <-function(n){return(rt(n,8)/sqrt(8/(8-2)))}

superposition_3graphiques(condition_stationnarite(rnorm),condition_stationnarite(runif_normalisee),condition_stationnarite(rt_8_normalisee))
ggsave("Condition de stricte stationnarité",path="./Graphiques_pour_Latex")



