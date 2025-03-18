# Chapitre IV : Intégration numérique

Ce chapitre porte sur les méthodes numériques pour l'approximation de l'intégrale d'une fonction.

---

## Position du problème

### Motivation

### Exemple de problème

Au cours de ce chapitre, nous appliquerons les différentes méthodes d'intégration à un même exemple : **La modélisation de la réfléctivité radar des gouttes de pluies**.

En 1948, Marshall et Palmer ont proposé un modèle du facteur de réflectivité des gouttes de pluies $Z$ (en $mm^6 m^{-3}$) pour les radars météorologiques :

$Z = \int_{0}^{D_{max}} N_0 e^{- \Lambda D} D^6 dD$

Les paramètres de ce modèle sont :

- $D_{max}$ la plus grande taille de goutte, que nous fixerons à $6 mm$.

- $\Lambda$ une constante empirique en $mm^{-1}$, pour lequel Marshall et Palmer proposent $4.1 R^{-0.21}$, avec $R$ le taux de pluie que nous fixerons à $5 mm.h^{-1}$.

- $N_0$ une constante empirique en $m^{-3} mm^{-1}$ nommée "paramètre de forme", pour lequel Marshall et Palmer proposent $N_0 = 8000$.

Ce modèle est encore aujourd'hui utilisé pour l'interprétation des mesures des radars météorologiques, dans le but d'estimer les précipitations aux sol à partir des réflectivités mesurées.

Nous essayerons donc de calculer l'intégrale entre $x=0$ et $x=6$ de la fonction $f(x) = N_0 e^{- \Lambda x} x^6$, dont la valeur est d'environ $3146.31 mm^6/m^3$.

![Graphique de f](img/Chap4_exemple_fonction.png)

Sous Python on utilisera la bibliothèque Numpy :

~~~
import numpy as np
~~~

Puis, on définira les variables globales suivantes :

~~~
N0 = 8000 #Paramètre de forme (m3/mm)
R = 5 #Taux de pluie (mm/h)
~~~

La fonction $f$ sera définie comme :

~~~
def f(D):
    
    return N0*np.exp(-(4.1*R**-0.21)*D)*D**6
~~~
