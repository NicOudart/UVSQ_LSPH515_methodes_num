# Chapitre III : Interpolation polynomiale

Ce chapitre porte sur les méthodes numériques pour l'approximation de valeures inconnues d'une fonction à partir des valeurs connues, par interpolation avec un polynôme.

---

## Position du problème

### Motivation

Soit une fonction $f(x)$ connue en seulement $n+1$ points, appelés **points** ou **noeuds d'interpolation** $(x_i,f(x_i))$ avec $i=0,1,2,...,n$ de l'intervalle $[a,b]$.
Peut-on approcher $f(x)$ pour tout $x$ de $[a,b]$ par une fonction ?

Il existe une infinité de fonctions d'interpolation, mais le plus simple est d'approcher la fonction par un **polynôme de degré suffisamment élevé** pour que sa courbe passe par les point d'interpolation.

Un polynôme $p$ de degré inférieur ou égal à $n$ s'exprime dans la **base canonique** ${1,x,x^2,...,x^n}$ de la manière suivante :

$p(x) = \displaystyle\sum_{k=0}^{n} a_k x^k$ où les $a_k$ sont les coefficients du polynôme.

**L'interpolation polynomiale consiste donc à déterminer les coefficients $a_k$ tels que $p(x_i) = f(x_i)$ pour $i=0,1,2,...,n$.**

**NB :** L'interpolation polynomiale pourra sont utiles pour les méthodes numériques de calcul d'intégrales et de dérivées.

**Attention !** L'interpolation polynomiale et l'approximation polynomiales sont des approches différentes.
L'approximation polynomiale de données bruitées cherche un polynôme de degré inférieur au nombre de données qui ne passe pas nécessairement par tous les points connus.

### Existence et unicité d'un polynôme d'interpolation

|Théorème d'Evariste Galois|
|:-|
|Un polynôme de degré $n$ a **au plus** $n$ racines qui peuvent être réelles ou complexes conjuguées.|

D'où le corollaire :

|Corollaire|
|:-|
|On ne peut faire passer par $n+1$ points distincts **qu'un seul** polynôme de degré n.|

**Toutes les méthodes présentées dans la suite de ce chapitre doivent donc aboutir au même polynôme**.

### Exemple de problème

Au cours de ce chapitre, nous appliquerons les différentes méthodes numériques d'interpolation polynomiale à un même exemple : **l'estimation de la durée du jour à l'UFR des sciences de l'UVSQ**.

La durée du jour (temps entre le lever et le coucher du soleil) en heures peut être approximée par la fonction $f$ suivante :

$f(x) = \frac{48}{2 \pi} \arccos(\tan(\lambda) \tan(\arcsin(\sin(\frac{2 \pi x}{365}) \sin(\delta))))$

avec $x$ le jour depuis l'équinoxe de printemps, $\lambda$ la latitude du lieu, et $\delta$ la latitude des tropiques.

On sait que $\delta \approx 23.438403°$.
On prendra ici l'exemple de l'UFR des sciences de l'UVSQ, dont la latitude est $\lambda \approx 48.81094°$.

On peut voir que cette formule implique 6 appels à des fonctions trigonométriques, ce qui peut rendre non négligeable le temps nécessaire pour évaluer $f$ en un grand nombre de points.
D'où l'intérêt de n'évaluer la fonction qu'en un nombre limité de points, et d'utiliser l'interpolation polynomiale pour déterminer d'autres valeurs de $f(x)$.

Nous choisirons ici d'évaluer la fonction pour les 10 valeurs de $x$ suivantes :

|x = jours depuis l'équinoxe de printemps|30   |60  |90  |120 |150  |180  |240  |270  |300  |330  |
|:---------------------------------------|:---:|:--:|:--:|:--:|:---:|:---:|:---:|:---:|:---:|:---:|
|f(x) = durée du jour en heures          |10.24|8.73|8.04|8.63|10.09|11.84|15.16|15.95|15.47|14.06|

Et nous essayerons d'estimer la valeur de $f$ pour $x = 210$ (i.e. la durée du jour en heures pour le 210ème jour depuis l'équinoxe de printemps) par interpolation polynomiale.

![Graphique de f](img/Chap3_exemple_fonction.gif)

Sous Python on utilisera la bibliothèque Numpy :

~~~
import numpy as np
~~~

Puis, on définira les variables globales suivantes :

~~~
l = 48.81094*np.pi/180 #Latitude de l'UFR des sciences
a = 23.438403*np.pi/180 #Latitude des tropiques
~~~

La fonction $f$ sera définie comme :

~~~
def f(d):
    
    return 48/(2*np.pi)*np.arccos(np.tan(l)*np.tan(np.arcsin(np.sin(a)*np.sin(d*2*np.pi/365.25))))
~~~

On calculera alors les 10 valeurs "connues" de la fonction de la manière suivante :

~~~
x = np.array([30,60,90,120,150,180,240,270,300,330],dtype='float')
y = f(x)
~~~

(On a ici définit avec Numpy un vecteur de 10 valeurs $x$, auquel on a appliqué $f$. Le vecteur résultant est stocké dans $y$).

## Matrices de Vandermonde

Une 1ère approche pour obtenir un polynôme d'interpolation est la suivante : déterminer les coefficients $a_i$ du polynôme en résolvant les $(n+1)$ équations de collocation $p(x_i) = f(x_i)$ pour $i=0,1,2,...,n$.

Ceci revient à résoudre le système linéaire de $(n+1)$ équations à $(n+1)$ inconnues :

$\begin{cases}
a_0 + a_1 x_0 + a_2 x_0^2 + ... + a_n x_0^n = f(x_0)\\
a_0 + a_1 x_1 + a_2 x_1^2 + ... + a_n x_1^n = f(x_1)\\
...\\
a_0 + a_1 x_n + a_2 x_n^2 + ... + a_n x_n^n = f(x_n)
\end{cases}$

Ce système admet une unique solution si les $x_i$ sont distincts 2 à 2.

Il peut s'écrire sous la forme matricielle suivante :

$\begin{pmatrix}
  1 & x_0 & x_0^2 & \cdots & x_0^n \\
  1 & x_1 & x_1^2 & \cdots & x_1^n \\
  \vdots  & \vdots  & \ddots & \vdots  \\
  1 & x_n & x_n^2 &\cdots & x_n^n 
 \end{pmatrix}
 \begin{pmatrix}
  a_0\\
  a_1\\
  \vdots\\
  a_n 
 \end{pmatrix}
 =
 \begin{pmatrix}
  f(x_0)\\
  f(x_1)\\
  \vdots\\
  f(x_n) 
 \end{pmatrix}$
 
On reconnait ici une **matrice de Vandermonde**.

La matrice de Vandermonde est inversible si et seulement si les $x_i$ sont distincts 2 à 2.

Dans la pratique, l'inversion de la matrice de Vandermonde **ne conduit pas à une solution satisfaisante** :

- Le nombre d'équations / d'inconnues croit avec le nombre de points d'interpolations, augmentant le nombre d'opérations de l'ordre de $2 n^3 / 3$.

- Ce type de système est souvent mal conditionné, ce qui rend la solution numérique très sensible aux erreurs d'arrondi.

C'est pourquoi dans la suite, on va préférer des techniques exprimant le polynôme dans **une autre base que la base canonique**.

## Polynômes de Lagrange

### L'algorithme

|Base de Lagrange|
|:-|
|Soient des $x_i$ (avec $i=0,1,2,...,n$) 2 à 2 distincts.|
|On appelle base de Lagrange relative aux points x_i les polynômes :|
|$L_i(x) = \displaystyle\prod_{j=0 , j \neq i}^{n} \frac{(x-x_j)}{x_i-x_j}$|
|soit $L_i(x) = \frac{(x-x_0)(x-x_1)...(x-x_n)}{(x_i-x_0)(x_i-x_1)...(x_i-x_n)}$|
|$L_i$ vérifie $L_i(x_i)=1$ et $L_i(x_j)=0$ si $j \neq i$.|

La famille des $(L_i(x))$ forme une base de l'ensemble des polynômes, et le polynôme qui interpoles les valeurs de $f(x_i)$ aux points $x_i$ s'écrit :

$p(x) = \displaystyle\sum_{i=0}^{n} f(x_i) L_i(x) = f(x_0) L_0(x) + f(x_1) L_1(x) + ... + f(x_n) L_n(x)$

### Exemple

## Polynômes de Newton

## Erreur d'interpolation

## Interpolation aux noeuds de Chebychev

## Interpolation par morceaux

### Interpolation affine

### Interpolation par fonctions splines
