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

Sous Python, on utilisera la bibliothèque Numpy pour définir / manipuler des matrices :

~~~
import numpy as np
~~~

### Solution, rang et déterminant

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

Dans le cas où la matrice $A$ est carrée, on a même le théorème suivant :

|Théorème|
|:-|
|Lorsque la matrice $A$ est carrée de dimension $n \times n$,|
|avec $n$ la taille du vecteur $x$,|
|le système linéaire $A x = b$ admet une **unique solution** si et seulement si :|
|le **déterminant de $A$** noté $det(A)$ est non nul.|
|($A$ est alors inversible, et on note $A^{-1}$ son inverse)|

On dit alors que le **système est de Cramer** et on peut écrire :

$x = A^{-1} b$

Le **système homogène** $A x = 0$ admet toujours le vecteur nul comme solution. 
- Si $det(A) \neq 0$ c'est l'unique solution. 
- Sinon il y en a une infinité.

Dans le cas d'un **système non homogène** ($b \neq 0$), si $det(A) = 0$ :
- Soit $rang(A) = rang([A \ mid b])$ alors il y a une infinité de solutions.
- Soit $rang(A) \neq rang([A \mid b])$ alors il n'y a pas de solution.

### Conditionnement

Même quand un système linéaire admet une solution unique, cette solution peut ne pas être "stable".

Un système est dit "**mal conditionné**" si la solution est extrêmement sensible aux perturbations des coefficients $A + \Delta A$ et des seconds membres $b + \Delta b$.

Un déterminant petit est souvent indicateur d'un mauvais conditionnement.

Pour quantifier la sensibilité de l'erreur relative sur la solution du système linéaire $A x = b$ aux variation de $A$ et $b$, on peut estimer le **conditionnement de la matrice A** :

$\kappa(A) = \|A\| \|A^{-1}\|$

avec une norme matricielle à définir.

L'erreur relative sur la solution est inférieure à l'erreur relative sur les données multiplisée par $\kappa(A)$ :
- Si $\kappa(A)$ est petit (de l'ordre de l'unité) on dit que le conditionnement est bon.
- Si $\kappa(A) >> 1$ le système est dit mal conditionné.

Le calcul du conditionnement dépend du choix de la norme :

- La **norme infinie** : $\|A\|_{\infty} = max_{1 \leq i \leq n} \displaystyle\sum_{j=1}^{n} |a_{i,j}|$

- La **norme 1** : $\|A\|_1 = max_{1 \leq j \leq n} \displaystyle\sum_{i=1}^{n} |a_{i,j}|$

- La **norme 2** : $\|A\|_2 = \|A^T\|_2 = \sqrt{\rho(A A^T)} = \sqrt{\rho(A^T A)}

$rho(A)$ est le **rayon spectral** de $A$, que l'on définit comme : 
$rho(A) = max_{1 \leq i \leq n} |\lambda_i|
avec $\lambda_i$ les **valeurs propres** de $A$.

On utilisera surtout le conditionnement de $A$ au sens de la norme 1 et de la norme 2 :

$\kappa_1(A) = \|A\|_1 \|A^{-1}\|_1$

$\kappa_2(A) = \|A\|_2 \|A^{-1}\|_2$

Pour une matrice carrée $A$ d'ordre $n$ inversible, le conditionnement vérifie les propriétés suivantes :

- $\kappa(A) /geq 1$

- $\forall \alpha \in \mathbb{R}$, $\kappa(\alpha A) = \kappa(A)$

- $\kappa(A) = \kappa(A^{-1})$

- Si on note $\sigma_{min}^2$ et $sigma_{max}^2$ la plus petite et la plus grande valeur propre de $A A^T$ : $\kappa_2(A) = \frac{\sigma_{max}}{\sigma_{min}}$

- Si $A$ est une matrice **réelle symétrique** ($A = A^T$), si $lambda_{min}$ et $lambda_{max}$ sont la plus petite et la plus grande valeur propre de $A$ en valeur absolue, on a : $\kappa_2(A) = \mid \frac{\lambda_{max}}{\lambda_{min}} \mid$

- Si $A$ est une matrice **orthogonale** ($A A^T = A^T A = I$) alors $\kappa_2(A) = 1$

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

Admettons que le récepteur GPS se trouve aux coordonnées ECEF $(x_r,y_r,z_r) = (4205,158,4777)$, correspondant approximativement à la position de l'UFR des Sciences de l'UVSQ.

Si la position des 4 satellites est :

$(x_{s1},y_{s1},z_{s1}) = (14000,4000,25000)$
$(x_{s2},y_{s2},z_{s2}) = (24000,6000,15000)$
$(x_{s3},y_{s3},z_{s3}) = (9000,-14000,21000)$
$(x_{s4},y_{s4},z_{s4}) = (10000,16000,19000)$

Alors le temps de retard associé à chaque satellite sera approximativement :

$t_1 = 0.0759878 s$
$t_2 = 0.0767739 s$
$t_3 = 0.0735321 s$
$t_4 = 0.0735485 s$

(Il est à noter que les positions des satellites ont été choisies pour être réalistes des satellites GPS).

Le système devient alors :

$\begin{pmatrix}
  10000 & 2000 & -10000 \\
  -5000 & -18000 & -4000 \\
  -4000 & 12000 & -6000
 \end{pmatrix}
 \begin{pmatrix}
  x_r\\
  y_r\\
  z_r 
 \end{pmatrix}
 =
 \begin{pmatrix}
  -5404000\\
  -42977000\\
  -43586000
 \end{pmatrix}$ 
 
C'est ce système d'équations linéaires que nous chercherons à résoudre pour essayer de retrouver la position $(x_r,y_r,z_r)$ du récepteur.

## La règle de Cramer

### Théorème

La **règle de Cramer** est un théorème d'algèbre linéaire qui donne la solution d'un système de Cramer :

|Théorème de la règle de Cramer|
|:-|
|Le système de Cramer $A x = b$ avec|
|$A$ matrice carrée de taille $n \times n$|
|$x$ vecteur de taille $n$|
|$b$ vecteur de taille $n$|
|admet une solution **si et seulement si** $A$ est inversible.|
|Cette solution est donnée par :|
|$x_i = \frac{det(A_i)}{det(A)}$ pour $i=1,...,n$|
|où $A_i$ est la matrice carrée formée en remplaçant la i-ème colonne de $A$ par le vecteur $b$.|

Lorsque le système n'est pas de Cramer (donc si $det(A)=0$) :
- Si le déterminant d'une des racines $A_i$ est nul alors le système n'a pas de solution.
- La réciproque est fausse : il peut arriver qu'un système n'ait pas de solution alors que tous les $det(A_i)$ sont non-nuls.

Cette méthode est très couteuse en nombre d'opérations et devient donc inapplicable à de grands systèmes.

### Algorithme

### Exemple

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