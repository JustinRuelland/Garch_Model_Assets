###--------------- Puissance du test --------------------
# Pas le calcul de la puissance du test, 
# mais vérification pour des cas précis que le test rejette bien 

source(file="./simulation_series.R")
source(file= "./QML_Variance.R",local=TRUE)
source(file= "./prevision.R",local=TRUE)

# Changement de GARCH
n = 5000 #nombres de jours
theta2 = c(0.0001,0.12,0.85)
theta2 = c(0.0001,0.12,0.85)


N = 400 #nombre de alphas
M = 10**4 #nombre de betas 

a = seq(0.01,0.2,length.out=N)
b = seq(0.6,0.98,length.out=M)
y = array(0,dim=c(M,N))

mes_alpha = c()
mes_beta = c()
mes_p_val = c()

for(alpha in a){
  print(alpha)
  for(beta in b){
    theta2 = c(0.0001,alpha,beta)
    simulation_rendements_avec_changement_GARCH(5000,theta1,theta2)
    
    p_val = func_backtest(rendements,-1.96,1.96,8*n%/%10)$p.value
    if(p_val>0){
      mes_alpha = c(mes_alpha,alpha)
      mes_beta = c(mes_beta,beta)
      mes_p_val = c(mes_p_val,p_val)
    }
  }
}

M=as.data.frame(y)

p = ggplot(data = y)+geom_bin2d()
p
