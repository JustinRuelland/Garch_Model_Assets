library(ggplot2)
library(signal)
library(dplyr)
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

superposition_3graphiques_condi_statio<-function(df1,df2,df3){
  df = bind_rows("Loi normale" = df1, "Uniforme normalisée" = df2, "Student (8DDL) normalisée" = df3, .id = "Loi_des_eta")
  
  Graph_condi_satio = ggplot(df) + 
    geom_line(aes(x=a,y=y,color = Loi_des_eta)) +
    xlab("Alpha") +
    ylab ("Bêta")
  
  plot(Graph_condi_satio)
  ggsave("Condition de stationnarité.png", width = 10, height = 5,path="./Graphiques_pour_Latex/")
}


# Exemple d'utilisation de la deuxième fonction (qui repose sur la première)

# superposition_3graphiques_condi_statio(condition_stationnarite(rnorm),condition_stationnarite(runif_normalisee),condition_stationnarite(rt_8_normalisee))



