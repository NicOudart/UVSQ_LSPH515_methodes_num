# Chapitre IV : Intégration numérique

Ce chapitre porte sur les méthodes numériques pour la résolution d'un système linéaire d'équations.

---

## Position du problème

### Motivation

Notre but est de résoudre un **système linéaire** : un ensemble d'équations portant sur les mêmes inconnues.

Un système de $m$ équations linéaires à $n$ inconnues peut s'écrire sous la forme suivante :

$\begin{cases}
a_{1,1} x_1 + a_{1,2} x_2 + ... + a_{1,n} x_n = b_1\\
a_{2,1} x_1 + a_{2,2} x_2 + ... + a_{2,n} x_n = b_2\\
...\\
a_{m,1} x_1 + a_{m,2} x_2 + ... + a_{m,n} x_n = b_m
\end{cases}$

avec $x_1,x_2,...,x_n$ les **inconnues** et $a_{i,j}$ les **coefficients** du système.

On peut également écrire ce système sous forme matricielle :

$A x = b$

avec $A$ une matrice de coefficients réels de taille $m \times n$, $x$ est un vecteur de taille $n$ contenant les variables réelles recherchées, et $b$ est un vecteur contenant $m$ réels.

$A = 
 \begin{pmatrix}
  a_{1,1} & a_{1,2} & \cdots & a_{1,n} \\
  a_{2,1} & a_{2,2} & \cdots & a_{2,n} \\
  \vdots  & \vdots  & \ddots & \vdots  \\
  a_{m,1} & a_{m,2} &\cdots & a_{m,n} 
 \end{pmatrix}$
 
$x =
 \begin{pmatrix}
  x_1\\
  x_2\\
  \vdots\\
  x_n 
 \end{pmatrix}$
 
$b =
 \begin{pmatrix}
  b_1\\
  b_2\\
  \vdots\\
  b_n 
 \end{pmatrix}$
 
Si $m>n$, on dit le système **sur-déterminé**.
Si $m<n$, on dit le système **sous-déterminé**.

### Solution, rang, déterminant, conditionnement

Face à un système linéaire, il y a 3 cas possibles :

- Le système n'a pas de solution.
- Le système a une infinité de solution.
- Le système a une solution unique.

On peut savoir dans quel cas on se trouve avec le **rang de la matrice $A$**.

|Définition|
|:-|
|Le **rang** d'une matrice $A$ est le nombre de vecteurs lignes ou colonnes linéairement indépendants.|

Si $A$ est de dimensions $m \times n$, alors $rang(A) \leq min(m,n)$.

|Théorème de Rouché-Fontené|
|:-|
|Le système linéaire $A x = b$ avec|
|$A$ une matrice de taille $m \times n$,|
|$x$ un vecteur de taille $n$,|
|et $b$ un vecteur de taille $m$,|
|admet une solution **si et seulement si** :|
|$rang(A) = rang([A \mid b])$|
|Si de plus, $rang(A) = n$, alors le système admet une **unique solution**.|
|Sinon, le système admet une infinité de solutions.|

### Règle de Cramer

### Exemple de problème

Au cours de ce chapitre, nous appliquerons les différentes méthodes d'intégration à un même exemple : **Le positionnement par satellites GPS**.

Nous exprimerons ici les positions en km, avec des coordonnées dans le repère cartésien ECEF (Earth-Centered Earth-Fixed), ayant pour origine le centre de la Terre.

* Soit un récepteur au sol dont on veut connaitre la position $(x_r,y_r,z_r)$ dans ce repère cartésien.

* Soient 4 satellites de la constellation GPS, dont la position est connue dans ce même repère : $(x_{s1},y_{s1},z_{s1})$, $(x_{s2},y_{s2},z_{s2})$, $(x_{s3},y_{s3},z_{s3})$ and $(x_{s4},y_{s4},z_{s4})$.

* Chaque satellite émet un signal, qui est reçu avec un certain temps de retard par le récepteur. Ces temps de retard $(t_1,t_2,t_3,t_4)$ sont mesurés par le récepteur.

* On admet que les signaux émis par chaque satellite se déplacent à vitesse constante jusqu'au récepteur : $c = 3.10^5 km/s$.

La distance euclidienne entre chaque satellite et le récepteur doit être égale au temps de retard du signal multiplié par sa vitesse.
On en déduit facilement que les différentes variables sont liées par le système de 4 équations suivant :

$\begin{cases}
(x_r-x_{s1})^2 + (y_r-y_{s1})^2 + (z_r-z_{s1})^2 = (c t_1)^2\\
(x_r-x_{s2})^2 + (y_r-y_{s2})^2 + (z_r-z_{s2})^2 = (c t_2)^2\\
(x_r-x_{s3})^2 + (y_r-y_{s3})^2 + (z_r-z_{s3})^2 = (c t_3)^2\\
(x_r-x_{s4})^2 + (y_r-y_{s4})^2 + (z_r-z_{s4})^2 = (c t_4)^2
\end{cases}$

Que l'on peut développer :

$\begin{cases}
x_r^2 - 2 x_{s1} x_r + x_{s1}^2 + y_r^2 - 2 y_{s1} y_r + y_{s1}^2 + z_r^2 - 2 z_{s1} z_r + z_{s1}^2 = (c t_1)^2\\
x_r^2 - 2 x_{s2} x_r + x_{s2}^2 + y_r^2 - 2 y_{s2} y_r + y_{s2}^2 + z_r^2 - 2 z_{s2} z_r + z_{s2}^2 = (c t_2)^2\\
x_r^2 - 2 x_{s3} x_r + x_{s3}^2 + y_r^2 - 2 y_{s3} y_r + y_{s3}^2 + z_r^2 - 2 z_{s3} z_r + z_{s3}^2 = (c t_3)^2\\
x_r^2 - 2 x_{s4} x_r + x_{s4}^2 + y_r^2 - 2 y_{s4} y_r + y_{s4}^2 + z_r^2 - 2 z_{s4} z_r + z_{s4}^2 = (c t_4)^2
\end{cases}$

En soustrayant la 1ère équation aux 3 autres, on réduit le système à 3 équations :

$\begin{cases}
x_{s2}^2 + y_{s2}^2 + z_{s2}^2 - 2 x_{s2} x_r - 2 y_{s2} y_r - 2 z_{s2} z_r - x_{s1}^2 - y_{s1}^2 - z_{s1}^2 + 2 x_{s1} x_r + 2 y_{s1} y_r + 2 z_{s1} z_r = (c t_2)^2-(c t_1)^2\\
x_{s3}^2 + y_{s3}^2 + z_{s3}^2 - 2 x_{s3} x_r - 2 y_{s3} y_r - 2 z_{s3} z_r - x_{s1}^2 - y_{s1}^2 - z_{s1}^2 + 2 x_{s1} x_r + 2 y_{s1} y_r + 2 z_{s1} z_r = (c t_3)^2-(c t_1)^2\\
x_{s4}^2 + y_{s4}^2 + z_{s4}^2 - 2 x_{s4} x_r - 2 y_{s4} y_r - 2 z_{s4} z_r - x_{s1}^2 - y_{s1}^2 - z_{s1}^2 + 2 x_{s1} x_r + 2 y_{s1} y_r + 2 z_{s1} z_r = (c t_4)^2-(c t_1)^2
\end{cases}$

Qui revient en regroupant les termes :

$\begin{cases}
x_{s2}^2 - x_{s1}^2 + y_{s2}^2 - y_{s1}^2 + z_{s2}^2 - z_{s1}^2 - 2 (x_{s2}-x_{s1}) x_r - 2 (y_{s2}-y_{s1}) y_r - 2 (z_{s2}-z_{s1}) z_r = (c t_2)^2-(c t_1)^2\\
x_{s3}^2 - x_{s1}^2 + y_{s3}^2 - y_{s1}^2 + z_{s3}^2 - z_{s1}^2 - 2 (x_{s3}-x_{s1}) x_r - 2 (y_{s3}-y_{s1}) y_r - 2 (z_{s3}-z_{s1}) z_r = (c t_3)^2-(c t_1)^2\\
x_{s4}^2 - x_{s1}^2 + y_{s4}^2 - y_{s1}^2 + z_{s4}^2 - z_{s1}^2 - 2 (x_{s4}-x_{s1}) x_r - 2 (y_{s4}-y_{s1}) y_r - 2 (z_{s4}-z_{s1}) z_r = (c t_4)^2-(c t_1)^2
\end{cases}$

Et en réarrangeant on obtient finalement le système linéaire :

$\begin{cases}
(x_{s2}-x_{s1}) x_r + (y_{s2}-y_{s1}) y_r + (z_{s2}-z_{s1}) z_r = \frac{1}{2} (x_{s2}^2 - x_{s1}^2 + y_{s2}^2 - y_{s1}^2 + z_{s2}^2 - z_{s1}^2 - (c t_2)^2 + (c t_1)^2)\\
(x_{s3}-x_{s1}) x_r + (y_{s3}-y_{s1}) y_r + (z_{s3}-z_{s1}) z_r = \frac{1}{2} (x_{s3}^2 - x_{s1}^2 + y_{s3}^2 - y_{s1}^2 + z_{s3}^2 - z_{s1}^2 - (c t_3)^2 + (c t_1)^2)\\
(x_{s4}-x_{s1}) x_r + (y_{s4}-y_{s1}) y_r + (z_{s4}-z_{s1}) z_r = \frac{1}{2} (x_{s4}^2 - x_{s1}^2 + y_{s4}^2 - y_{s1}^2 + z_{s4}^2 - z_{s1}^2 - (c t_4)^2 + (c t_1)^2)
\end{cases}$

avec 3 équations et 3 inconnues $(x_r,y_r,z_r)$.

On peut écrire ce système sous la forme matricielle $A x = b$ :

$\begin{pmatrix}
  x_{s2}-x_{s1} & y_{s2}-y_{s1} & z_{s2}-z_{s1} \\
  x_{s3}-x_{s1} & y_{s3}-y_{s1} & z_{s3}-z_{s1} \\
  x_{s4}-x_{s1} & y_{s4}-y_{s1} & z_{s4}-z_{s1}
 \end{pmatrix}
 \begin{pmatrix}
  x_r\\
  y_r\\
  z_r 
 \end{pmatrix}
 =
 \begin{pmatrix}
  \frac{1}{2} (x_{s2}^2 - x_{s1}^2 + y_{s2}^2 - y_{s1}^2 + z_{s2}^2 - z_{s1}^2 - (c t_2)^2 + (c t_1)^2)\\
  \frac{1}{2} (x_{s3}^2 - x_{s1}^2 + y_{s3}^2 - y_{s1}^2 + z_{s3}^2 - z_{s1}^2 - (c t_3)^2 + (c t_1)^2)\\
  \frac{1}{2} (x_{s4}^2 - x_{s1}^2 + y_{s4}^2 - y_{s1}^2 + z_{s4}^2 - z_{s1}^2 - (c t_4)^2 + (c t_1)^2)
 \end{pmatrix}$
 
C'est ce système d'équations linéaires que nous chercherons à résoudre pour déterminer la position $(x_r,y_r,z_r)$ du récepteur.

Pour vérifier nos résultats, nous considérerons que le récepteur GPS se trouve aux coordonnées ECEF $(x_r,y_r,z_r) = (4205,158,4777)$, correspondant approximativement à la position de l'UFR des Sciences de l'UVSQ.

## Méthodes directes d'élimination

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