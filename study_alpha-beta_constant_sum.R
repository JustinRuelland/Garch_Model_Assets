# In the work on the power of the test for GARCH change (puissance_test_chgtGARCH.R), 
# we found (in the case where the change in training set and test set corresponds to the date of GARCH change), 
# the existence of a straight line with slope -1 where all pairs (alpha,beta) for the second GARCH 
# were not rejected by the chi-square test.

# Hence, we focus here on the understanding of the model for two pairs of (alpha,beta) 
# such that alpha+beta is constant.
rm(list=ls())

#---------------- import -------------------
source(file = "./simulation_series.R", local=TRUE)
source(file = "./QML_Variance.R", local=TRUE)
source(file= "./prevision.R", local=TRUE)


#-------------- Studied pairs (alpha,beta) ---------------
# theta1 et theta2 font ici référence aux deux theta étudiés pour la deuxiè
n = 1000 #days
cut = 0.8 #cut du backtest ET du changement de GARCH

omega = 0.0001

theta0 = c(omega, 0.12, 0.85)

theta1 = theta0

alpha2 = 0.8 #must be included in [0,0,97] (0.97 = alpha1+beta1)
theta2 = c(omega,alpha2,theta1[2]+theta1[3]-alpha2)
#theta2 = c(omega,0.12,0.12)

#------------------- simulation ------------------------
set.seed(3)
etas = rnorm(n)

eps_1 = simulation_rendements_avec_changement_GARCH(n,theta_1 = theta0, theta_2 = theta1, cut = cut, etas = etas)
eps_2 = simulation_rendements_avec_changement_GARCH(n,theta_1 = theta0, theta_2 = theta2, cut = cut, etas = etas)


courbe_prix_arg_rendements(eps_1)
courbe_prix_arg_rendements(eps_2)

#--------------------- STUDY -------------------

pval_1 = func_backtest(eps_1,-1.96,1.96,empirical = FALSE, cut = 0.8)$p.val #Remarque : le QML est refait dans la fonction
pval_2 = func_backtest(eps_2,-1.96,1.96,empirical = FALSE, cut = 0.8)$p.val

pval_1
pval_2

n_cut = floor(cut*n)
theta1_hat = QML(eps_1[1:n_cut]**2)
theta2_hat = QML(eps_2[1:n_cut]**2)

sigma1_hat = simu_sigma2(eps_1**2, theta=theta1_hat)
sigma2_hat = simu_sigma2(eps_2**2, theta=theta2_hat)

# Affichages - que de n_cut à n -> test set du backtest
window = n_cut:820

df = cbind(Date = window,sigma1_hat[window],sigma2_hat[window])
df = as.data.frame(df)
colnames(df) = c("Date","sigma1_hat","sigma2_hat")

p1 = ggplot(data=df)+geom_line(aes(x = Date,y = sigma1_hat))+ylab("Carré du sigma estimé (série 1)")
p1

p2 = ggplot(data=df)+geom_line(aes(x=Date, y = sigma2_hat))+ylab("Carré du sigma estimé (série 2)")
p2

# Comparaison des moyennes
m1 = mean(sigma1_hat)
m2 = mean(sigma2_hat)
# print("blabl" + m1)
