# Garch_Model_Assets
## Les fonctions implémentées
- Fichier **data_preparation.R** :
  - *transform_csv* : fonction pour transformer le csv de Yahoo finance en dataframe avec les rendements et sans les variables qui ne nous intéresse pas
  - *plot_series_temp* : réalise les graphiques des prix, rendements et rendements au carré en fonction du temps
  - *autocorrelations* : graphiques autocorrélations des rendements et rendements au carré

- Fichier **condition_stationnarite.R** :
  - *condition_stationnarite* : à partir de la loi des etas (bruit blanc), renvoie un dataframe de couples (alpha,beta) qui indiquent la condition de stationnarité
  - *superposition_3graphiques_condi_statio* : à partir de 3 dataframes (comme ceux de *condition_stationnarite*), renvoie le graphs des trois courbes de la condition de stationnarité

- Fichier **simulation_series.R** :
  - *courbe_prix_arg_rendements* : à partir d'un dataframe de rendements, affiche la courbe des prix (avec le prix en t=0 normalisé à 1)
  - *simulation_rendements* : à partir d'un nombre n de simulations, d'un paramètre theta, d'une loi des éta (par défaut, rnorm) et d'un sigma carré initial (par défaut,  $\omega=\frac{1}{1-\alpha-\beta}$), simule n carrés de rendement pour un GARCH(1,1) et renvoie sous forme de liste
  - *simulation_rendements_avec_changement_GARCH* : à partir d'un nombre n de simulations, de deux paramatères theta d'une loi des éta (par défaut, rnorm) et d'un sigma carré initial (par défaut, $\omega=\frac{1}{1-\alpha-\beta}$ avec les valeurs du premier theta), renvoie sous forme de liste la simulation d'un GARCH(1,1) avec changement de paramètre au milieu
  - *simulation_GARCH22* : à partir d'un nombre n de simulations, d'un THETA ( $=(\omega,\alpha_1,\alpha_2,\beta_1,\beta_2)$ , d'une loi des éta (par défaut, rnorm), et d'un sigma caré initial (par défaut, $\omega=\frac{1}{1-\alpha_1-\alpha_2-\beta_1-\beta_2}$ ), renvoie une liste d'une série GARCH(2,2) de n termes

- Fichier **QML_Variance.R**:
  - *simu_eps*: fonction pour simuler une suite de rendements (notés epsilon) à partir de la forume du modèle Garch. La longueur de la suite, les paramètres du Garch et les valeurs initiales sont laissées en paramètre.
  - *simu_sigma2*: fonction pour simuler la volalité sigma2 à partir des rendements au carré (*eps2*) observés et d'un paramètre theta. En pratique, on aura seulement les rendements au carré et on évaluera par QML un $\hat{\theta}$ pour ensuite simuler des $\hat{\sigma^{2}}$.
  - *QML*: fonction pour estimer par QML les paramètres *theta_hat* du modèle à partir d'observations *eps2*.
  - *var_asymp*: fonction pour estimer la matrice de variance asymptotique de la suite $\sqrt{n}(\hat{\theta}-\theta_{0}$, $\theta_{0}$ étant le vrai paramètre.
