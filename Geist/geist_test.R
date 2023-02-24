source(file = "./data_preparation.R",local= TRUE)
data = read.csv("./CAC40_15_19.csv") 
#data = read.csv("./^GDAXI.csv")

data_disc = transform_csv_with_discount(data,0.01/365)
plot_series_temp(data_disc)

data_no_discount = transform_csv(data)
plot_series_temp(data_no_discount)


# Comparaison rapide
source(file= "./QML_Variance.R",local=TRUE)
source(file= "./prevision.R",local=TRUE)

rendements_disc = c(data_disc$rendement) # si on met c() ou pas, ça change la taille du vecteur -> comprends pas
n = length(rendements_disc)
p_val_disc = func_backtest(rendements_disc,-1.96,1.96,8*n%/%10,FALSE)$p.value

rendements_no_disc = c(data_no_discount$rendement)
p_val_no_disc = func_backtest(rendements_no_disc,-1.0,1.0,8*n%/%10,FALSE)$p.value


rendements_disc[200]
rendements_no_disc[200]

print(p_val_disc == p_val_no_disc)

# Pb pour DAX : la fonction d'optimisation bug car série trop longue
# CAC40 : p_valeur nickel pour
