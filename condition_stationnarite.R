
#---------------------- Condition de stationnarité -------------------------------
### But : Représentations des (alpha,beta) vérifiant la condition de stationnarité

### Utilisation de la fonction :
# Mettre en paramètres la fonction de simulation et ses paramètres
# Par exemple "f = rnorm", et les paramètres 0 et 1

cond_statio <- function(f,param1,param2){
  N = 400 #nombre de alphas
  M = 10**4 #nombre de betas (ne pas dépasser 10**4 à priori)
  nsimu= 10**5 #nombre de simulations
  
  a = seq(0,2,length.out=N)
  b = seq(0,1,length.out=M)
  y = array(0,dim=c(M,N)) #résultat final
  
  for(i in range(nsimu)){
    #on évalue tous les alphas et betas sur une simulation
    xx = array(f(N*M,param1,param2), dim=c(N,M)) 
    y = y + log(t(a*xx**2)+b)/nsimu #on ajoute à la moyenne le résultat
  }
  
  #0 si espérance négative, 1 si positive
  y[y>0]=0
  y[y<0]=1
  y = apply(y,2,FUN=mean) # sur les betas
  
  df = cbind(a,y)
  df = as.data.frame(df)
  Graph_condi_satio = ggplot(data = df) + geom_line(aes(x = a,y = y))+scale_x_continuous("Alpha") + scale_y_continuous("Beta")
  plot(Graph_condi_satio)
  
}
