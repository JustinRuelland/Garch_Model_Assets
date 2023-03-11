# In the work on the power of the test for GARCH change (puissance_test_chgtGARCH.R), 
# we found (in the case where the change in training set and test set corresponds to the date of GARCH change), 
# the existence of a straight line with slope -1 where all pairs (alpha,beta) for the second GARCH 
# were not rejected by the chi-square test.

# Hence, we focus here on the understanding of the model for two pairs of (alpha,beta) 
# such that alpha+beta is constant.
rm(list=ls())

#---------------- import -------------------
source(file = "./simulation_series.R",local=TRUE)
#source(file = )


#-------------- Studied pairs (alpha,beta) ---------------
n = 1000 #days

omega = 0.0001

theta0 = c(omega,0.12,0.85)

theta1 = theta0

alpha2 = 0.8 #must be included in [0,0,97] (0.97 = alpha1+beta1)
theta2 = c(omega,alpha2,theta1[2]+theta1[3]-alpha2)

#------------------- simulation ------------------------
set.seed(2)
eps_1 = simulation_rendements(n,theta1)
eps_2 = simulation_rendements(n,theta2)

eps_1 = simulation_rendements_avec_changement_GARCH(n,theta_1 = theta0, theta_2 = theta1, cut = 0.8)

courbe_prix_arg_rendements(eps_1)
courbe_prix_arg_rendements(eps_2)

#--------------------- STUDY -------------------

