# Garch_Model_Assets
## Les fonctions implémentées
- Fichier **data_preparation.R** :
  - *transform_csv* : fonction pour transformer le csv de Yahoo finance en dataframe avec les rendements et sans les variables qui ne nous intéresse pas
  - *plot_series_temp* : réalise les graphiques des prix, rendements et rendements au carré en fonction du temps
  - *autocorrelations* : graphiques autocorrélations des rendements et rendements au carré
  - *transform_csv_with_discount*(data,r) : transform_csv mais qui discount les prix au taux quotidient r
  
- Fichier **white_noise_laws.R** : On implémente des lois normalisées pour la loi des eta, pour qu'on puisse mettre autre chose que "rnorm(n)". Toutes les fonctions renvoient un vecteur de taille n.
  - runif_normalisee(n) : loi uniforme normalisée
  - rt_8_normalisee(n) : loi de Student à 8 degrés de liberté normalisée
  - normalised_student(n) : loi de Student à 5 degrés de liberté normalisée (le degré ici a vocation à être modifiée si nécessaire)

- Fichier **condition_stationnarite.R** :
  - *condition_stationnarite* : à partir de la loi des etas (bruit blanc), renvoie un dataframe de couples (alpha,beta) qui indiquent la condition de stationnarité
  - *superposition_3graphiques_condi_statio* : à partir de 3 dataframes (comme ceux de *condition_stationnarite*), renvoie le graphs des trois courbes de la condition de stationnarité

- Fichier **simulation_series.R** :
  - *courbe_prix_arg_rendements* : à partir d'un dataframe de rendements, affiche la courbe des prix (avec le prix en t=0 normalisé à 1)
  - *simulation_rendements* : à partir d'un nombre n de simulations, d'un paramètre theta, d'une loi des éta (par défaut, rnorm) et d'un sigma carré initial (par défaut,  $\omega=\frac{1}{1-\alpha-\beta}$), simule n carrés de rendement pour un GARCH(1,1) et renvoie sous forme de liste
  - *simulation_rendements_avec_changement_GARCH* $(n,\theta_1,\theta_2, cut = 0.5, sigma\\_init = \sqrt{\frac{\theta_1[1]}{1-\theta_1[2]-\theta_1[3]}}), etas = rnorm(n))$ : à partir d'un nombre n de simulations, de deux paramatères theta, d'un vecteur de eta de taille n (par défaut, rnorm(n)) et d'un sigma carré initial (par défaut, $\omega=\frac{1}{1-\alpha-\beta}$ avec les valeurs du premier theta), renvoie sous forme de liste la simulation d'un GARCH(1,1) avec changement de paramètre au "cut" (cut $\in [0,1]$)
  - *simulation_GARCH22* : à partir d'un nombre n de simulations, d'un THETA ( $=(\omega,\alpha_1,\alpha_2,\beta_1,\beta_2)$ , d'une loi des éta (par défaut, rnorm), et d'un sigma caré initial (par défaut, $\omega=\frac{1}{1-\alpha_1-\alpha_2-\beta_1-\beta_2}$ ), renvoie une liste d'une série GARCH(2,2) de n termes

- Fichier **QML_Variance.R**:
  - *simu_eps*: fonction pour simuler une suite de rendements (notés $\epsilon$ dans les notes et *eps* dans le code) à partir de la formule du modèle Garch. La longueur de la suite, les paramètres du Garch et les valeurs initiales sont laissés en paramètre.
  - *simu_sigma2*: fonction pour simuler la volalité (notes: $\sigma^{2}$, code: *sigma2*) à partir des rendements au carré ( $\epsilon^{2}$, *eps2*) observés et d'un paramètre $\theta$. En pratique, on aura seulement les rendements au carré et on évaluera par QML un $\hat{\theta}$ pour ensuite simuler des $\hat{\sigma^{2}}$.
  - *QML*: fonction pour estimer par QML les paramètres $\hat{\theta}$ du modèle à partir d'observations de $\epsilon^{2}$.
  - *var_asymp*: fonction pour estimer la matrice de variance asymptotique de la suite $\sqrt{n}(\hat{\theta}-\theta_{0})$, $\theta_{0}$ étant le vrai paramètre.

- Fichier **prevision.R**:
  - *func_backtest*: d'abord cette fonction divise les données en deux parties en fonction d'un paramètre *cut*. Le modèle est estimé sur la première partie et un intervalle de confiance est construit pour chaque donnée de la second partie. La fonction renvoie les bornes des intervalles de confiance. De plus théoriquement, $\alpha$% des valeurs doivent être hors des intervalles de confiance, ceci est testé avec un test du khi2, la p-valeur du test est aussi renvoyée. $\alpha$ et les quantiles de la loi choisie pour modéliser $\eta$ sont laissés en paramètre.

- Fichier **puissance_test_chgtGARCH.R** :
  - *puissance_test_chgtGARCH* : simule (peut prendre facilement 20min quand n_path =25 pex, ou n=2000) une "carte bleue" (=représentation des zones de rejet de test dans le cadre d'une simulation GARCH(1,1) avec changement de GARCH à un "cut") en fonction de :
    - cut_chgt = 0.4 : part du temps simulé sous GARCH(1,1) theta_1
    - n_path = 10 : nombres de trajectoires simulées pour chaque couple $(\alpha,\beta)$
    - n = 1000: nombres de jours de simulation
    - Remarque : on utilise les quantiles théoriques de la loi normale (donc +-1.96) pour faire le test
    - backtest : spécifie la fonction de backtest qu'on utilise (du fichier prevision.R)

- Fichier **study_alpha-beta_constant_sum.R** :
  - *sigma2_hat_alpha_beta_constant* : le but de cette fonction est d'afficher quelques éléments autour d'une simulation pour un changement de GARCH avec $\alpha+\beta$ constant. En particulier, ça affiche :
    - quelques statistiques sur les trajectoires
    - courbe des prix des cas avec et sans changement de GARCH
    - graphe des simga2 réels sur la période $0.8 n$ à $n$ (contient aussi un changement GARCH hors de la droite $\alpha+\beta$ constant)
    - graphe des $\sigma^2$ et des $\hat{\sigma}^2 estimés sur la période $0.8 n$ à $n$
  - *test_puissance_changment_horizon_long* : pour montrer que y'a effectivement plus de rejets même sur la droite $\alpha+\beta$ constant lorsque la période est plus grande
