# Chapitre IV : Intégration numérique

Ce chapitre porte sur les méthodes numériques pour la résolution d'un système linéaire d'équations.

---

## Position du problème

## Méthodes directes d'élimination

### Motivation

### Solution, rang, déterminant, conditionnement

### Règle de Cramer

### Exemple de problème

Au cours de ce chapitre, nous appliquerons les différentes méthodes d'intégration à un même exemple : **Le positionnement par satellites GPS**.

* Soit un récepteur au sol dont on veut connaitre la position $(x_r,y_r,z_r)$ dans un repère cartésien.

* Soient 3 satellites de la constellation GPS, dont la position est connue dans ce même repère : $(x_{s1},y_{s1},z_{s1})$, $(x_{s2},y_{s2},z_{s2})$ et $(x_{s3},y_{s3},z_{s3})$.

* Chaque satellite émet un signal, qui est reçu avec un certain temps de retard par le récepteur. Ces temps de retard $(t_1,t_2,t_3)$ sont mesurés par le récepteur.

* On admet que les signaux émis par chaque satellites se déplacent à vitesse constante jusqu'au récepteur : $c = 3.10^8 m/s$.

La distance euclidienne entre chaque satellite et le récepteur doit être égale au temps de retard du signal multiplié par sa vitesse.
On en déduit facilement que les différentes variables sont liées par le système de 3 équations suivant :

$\begin{cases}
(x_r-x_{s1})^2 + (y_r-y_{s1})^2 + (z_r-z_{s1})^2 = (c t_1)^2\\
(x_r-x_{s2})^2 + (y_r-y_{s2})^2 + (z_r-z_{s2})^2 = (c t_2)^2\\
(x_r-x_{s3})^2 + (y_r-y_{s3})^2 + (z_r-z_{s3})^2 = (c t_3)^2
\end{cases}$

Que l'on peut développer :

$\begin{cases}
x_r^2 - 2 x_{s1} x_r + x_{s1}^2 + y_r^2 - 2 y_{s1} y_r + y_{s1}^2 + z_r^2 - 2 z_{s1} z_r + z_{s1}^2 = (c t_1)^2\\
x_r^2 - 2 x_{s2} x_r + x_{s2}^2 + y_r^2 - 2 y_{s2} y_r + y_{s2}^2 + z_r^2 - 2 z_{s2} z_r + z_{s2}^2 = (c t_2)^2\\
x_r^2 - 2 x_{s3} x_r + x_{s3}^2 + y_r^2 - 2 y_{s3} y_r + y_{s3}^2 + z_r^2 - 2 z_{s3} z_r + z_{s3}^2 = (c t_3)^2
\end{cases}$

Puis simplifier :

$\begin{cases}
x_{s1}^2 + y_{s1}^2 + z_{s1}^2 - 2 x_{s1} x_r - 2 y_{s1} y_r - 2 z_{s1} z_r = (c t_1)^2\\
x_{s2}^2 + y_{s2}^2 + z_{s2}^2 - 2 x_{s2} x_r - 2 y_{s2} y_r - 2 z_{s2} z_r = (c t_2)^2\\
x_{s3}^2 + y_{s3}^2 + z_{s3}^2 - 2 x_{s3} x_r - 2 y_{s3} y_r - 2 z_{s3} z_r = (c t_3)^2
\end{cases}$

On obtient finalement un système linéaire :

$\begin{cases}
x_{s1} x_r + y_{s1} y_r + z_{s1} z_r = \frac{1}{2} (x_{s1}^2 + y_{s1}^2 + z_{s1}^2 - (c t_1)^2)\\
x_{s2} x_r + y_{s2} y_r + z_{s2} z_r = \frac{1}{2} (x_{s2}^2 + y_{s2}^2 + z_{s2}^2 - (c t_2)^2)\\
x_{s3} x_r + y_{s3} y_r + z_{s3} z_r = \frac{1}{2} (x_{s3}^2 + y_{s3}^2 + z_{s3}^2 - (c t_3)^2)
\end{cases}$

avec 3 équations et 3 inconnues $(x_r,y_r,z_r)$.

On peut écrire ce système sous la forme matricielle $A x = b$ :

$\begin{pmatrix}
  x_{s1} & y_{s1} & z_{s1} \\
  x_{s2} & y_{s2} & z_{s2} \\
  x_{s3} & y_{s3} & z_{s3} 
 \end{pmatrix}
 \begin{pmatrix}
  x_r\\
  y_r\\
  z_r 
 \end{pmatrix}
 =
 \begin{pmatrix}
  \frac{1}{2} (x_{s1}^2 + y_{s1}^2 + z_{s1}^2 - (c t_1)^2)\\
  \frac{1}{2} (x_{s2}^2 + y_{s2}^2 + z_{s2}^2 - (c t_2)^2)\\
  \frac{1}{2} (x_{s3}^2 + y_{s3}^2 + z_{s3}^2 - (c t_3)^2)
 \end{pmatrix}$

### Propriétés des systèmes linéaires

### Pivot de Gauss

### Elimination de Gauss-Jordan

## Méthodes directes de factorisation / décomposition

### Equivalence élimination - factorisation

### Décomposition LU et PLU

### Autres décompositions (QR et Cholesky)

## Méthodes itératives

### Principe et convergence

### Méthode de Jacobi

### Méthode de Gauss-Seidel

### Méthode de relaxation

### Méthodes du gradient

## Conclusion