# Garch_Model_Assets
## Les fonctions implémentées
- Fichier **data_preparation.R** :
  - *transform_csv* : fonction pour transformer le csv de Yahoo finançant en dataframe avec les rendements et sans les variables qui ne nous intéresse pas
  - *plot_series_temp* : réalise les graphiques en fonction du temps des prix, rendements et rendements au carré
  - *autocorrelations* : graphiques autocorrélations des rendements et rendements au carré

- *condition_stationnarite.R* :
  - *condition_stationnarite* : à partir de la loi des etas (bruit blanc), renvoie un dataframe de couples (alpha,beta) qui indiquent la condition de stationnarité
  - *superposition_3graphiques_condi_statio* : à partir de 3 dataframes (comme ceux de *condition_stationnarite*), renvoie le graphs des trois courbes de la condition de stationnarité

- *simulation_series.R* :
  - *courbe_prix_arg_rendements* : à partir d'un dataframe de rendements, affiche la courbe des prix (avec le prix en t=0 normalisé à 1)
  - *simulation_rendements* : Simulation de n carrés de rendement pour un GARCH(1,1) de paramètres theta et avec un sigma carré initial (par défaut, le sigma carré initial est donnée par $\omega=\frac{1}{1-\alpha-\beta}$)

