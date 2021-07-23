# -tresthor-

Le package tresthor développé par la DG Trésor est un outil permettant de créer des modèles économiques, de les résoudre pour effectuer des prévisions et de les analyser. Il est mis à disposition du public afin d'éclairer sur les méthodes utilisées pour les prévisions économiques des projets de loi de finances et programmes de stabilités. Il peut être également utilisé à des fins pédagogiques pour mieux appréhender les mécanismes économiques suivant un modèle macroéconomique construit par l'utilisateur.

* Les dépendances du package :
  + `dplyr` et `Deriv`(seront chargés avec le package)
  + `purrr`, `tidyr`, `stringr`, `ggplot2`
  + `scales` , `splitstackshape`, `assertthat`,`stats`,`Matrix`,`cointReg`,`gsubfn`
  + `Rcpp`, `RcppArmadillo`
* Le package fonctionne sur R en version 4.0.2 ou ultérieure, et est compatible avec Microsoft Open R.

* Pour faire fonctionner le solveur en mode Rcpp (entre autres), il peut être nécessaire d'installer Rtools sur Windows, et sur MacOS il faudra [XQuartz](https://www.xquartz.org) et soit XCode ou une version de [gfortran](https://github.com/fxcoudert/gfortran-for-macOS) pour avoir un compilateur C++ compatible.

# Documentation disponible

- Le [guide de l'utilisateur](https://www.tresor.economie.gouv.fr/Content/other/Opale/manuel_tresthor.html) présente en détail l'ensemble des fonctionnalités du package.
- Le [guide d'utilisation du modèle Opale](https://www.tresor.economie.gouv.fr/Content/other/Opale/manuel_opale_r.html) en R par tresthor.
- Le [guide appliqué de prévision](https://www.tresor.economie.gouv.fr/Content/other/Opale/application_tresthor_modele_uk.html) par tresthor avec un modèle simplifié sur le Royaume-Uni.</i>

# Opale : les ressources mises à disposition

Opale est le modèle macroéconomique de la DG Trésor utilisé pour les prévisions de croissance à horizon 2-3 ans pour les budgets des lois de finances et programmes de stabilité. Il s'agit d'un modèle composé d'environ 500 équations pour représenter l'économie française à partir d'équations économétriques et comptables.

La documentation complète sur la construction de ce modèle est disponible sur le [site de la DG Trésor](https://www.tresor.economie.gouv.fr/Articles/2017/05/19/la-maquette-de-prevision-opale-2017).

  * La version d'Opale exploitable par `tresthor` sur R est incluse dans le package, accessible avec le chemin d'accès suivant : `system.file("Opale","opale.txt",package = "tresthor")`. Il faut noter que quelques équations ont subi des changements depuis la publication du document de travail mentionné ci-dessus.

  * Nous mettons à disposition un jeu de données trimestrielles allant de 1980 à 2020, basé sur les comptes trimestriels de l'Insee publiés le 28 mai 2021 et avec quelques données supplémentaires utilisées dans le modèle. Elles sont accessibles directement dans le package `tresthor` installé par le chemin suivant : `system.file("Opale","donnees_opale.rds",package = "tresthor")`.

  * Les coefficients économétriques estimés et non inscrits directement dans le modèles sont disponibles dans le fichier suivant : `system.file("Opale","donnees_opale.rds",package = "tresthor")`. 

  * Un dictionnaire des variables utilisées pour faciliter la compréhension du modèle est accessible dans le package installé : `View(dictionnaire_variables_opale)`.



  Ce package a été réalisé par Anissa Saumtally, avec la collaboration de Charlotte Nudelmann et Niamh Dunne. Nous remercions Pierre Lissot et Antonin Aviat, ainsi que tous les testeurs et relecteurs qui ont aidé au développement de ce package.

  Ce package en logiciel libre est régi par une licence CeCILL-C. Le code source est disponible ici.

  # Installation du package

  1. Télécharger la [dernière version du package](https://framagit.org/DGTresor/opale-en-r-tresthor/-/blob/main/tresthor_1.0.0.tar.gz) sur ce dépot framagit (au 29/06/2021 : tresthor_1.0.0.tar.gz)
  2. Installer les dépendances du package dans R avec la commande `install.packages(c("Deriv","tidyverse","scales","splitstackshape","assertthat","Matrix","cointReg","Rcpp","RcppArmadillo","gsubfn"))`
  3. Installer le package téléchargé :

  * manuellement dans RStudio : en passant par Tools/Install packages, puis en choisissant l'option "Package Archive File (.zip, .tar.gz)" et enfin sélectionner le fichier .tar.gz téléchargé)

  * avec la ligne de commande `install.packages("chemin du fichier/tresthor_1.0.0.tar.gz", repos=NULL)`
   
  
  # Code source

  Le code source du package est disponible [ici sur le projet git tresthor_dev](https://framagit.org/DGTresor/tresthor_dev).

