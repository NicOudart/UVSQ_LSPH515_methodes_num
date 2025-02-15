# Chapitre II : Recherche de racines

Ce chapitre porte sur les méthodes numériques pour la recherche de racines d'une fonction.

---

## Position du problème

### Motivation

|Définition|
|:-|
|Soit f une fonction définie et continue de $\mathbb{R}$ dans $\mathbb{R}$.|
|Les **racines** ou **zéros** de cette fonction sont les valeurs $x$ qui vérifient l'équation $f(x)=0$|

Le principal intérêt des méthodes de recherche de racines est de pouvoir **résoudre des équations**.
En effet, trouver $x \in \mathbb{R}$ tel que $f(x)=c$ revient à chercher les racines de la fonction $f(x)-c$.

Dans certains cas, on peut trouver analytiquement les racines d'une fonction.  

Toute équation polynomiale de degré $n$ a exactement $n$ solutions dans $\mathbb{C}$ et au plus $n$ solutions dans $\mathbb{R}$.
Mais d'après la théorie d'Evariste Galois, **à partir du degré 5 il n'existe plus de formule générale de résolution**.

Si l'équation n'est pas polynomiale, sa résolution analytique est encore moins probable.

D'où l'intérêt d'utiliser des **méthodes de résolution numériques**, afin d'approcher les valeurs des racines.

Le principe est le suivant :

* Localiser grossièrement les racines en procédant à des évaluations graphiques.

* Construire une suite qui converge vers chaque racine.

### Existence et localisation des racines

Que l'approche soit analytique ou numérique, la 1ère étape consiste généralement à localiser les solution de l'équation.

Pour ce faire, on détermine les intervalles $[a,b]$ contenant une unique **racine**.
C'est ce que l'on appelle la **séparation des racines**.

Cette détermination se fait graphiquement et/ou analytiquement de la manière suivante :

* Etude des variations de la fonction $f$.

* Trouver les intervalles ne contenant qu'une seule racine en s'appuyant sur les théorèmes suivant.

|Théorème des valeurs intermédiaires|
|:-|
|Soit $f$ une fonction continue sur un intervalle $I=[a,b]$ de $\mathbb{R}$.|
|$f$ atteint toutes les valeurs intermédiaires entre $f(a)$ et $f(b)$. Autrement dit :|
|- Si $f(a) \leq f(b)$ alors pour tout $d \in [f(a),f(b)]$ il existe un $c \in [a,b]$ tel que $f(c)=d$.|
|- Si $f(a) \geq f(b)$ alors pour tout $d \in [f(b),f(a)]$ il existe un $c \in [a,b]$ tel que $f(c)=d$.|

D'où le corollaire :

|Théorème de Bolzano|
|:-|
|Soit une fonction continue $f:[a,b] \longrightarrow \mathbb{R}$| 
|Si $f(a)f(b)<0$ alors il existe au moins un $c \in ]a,b[$ tel que $f(c)=0$.|

Le théorème de Bolzano garantit l'existence d'au moins une racine **mais pas son unicité** dans $[a,b]$.

Pour assurer l'unicité d'une racine dans $[a,b]$, on essaiera d'appliquer le **théorème de la bijection** avec le **théorème des valeurs intermédiaires**.

|Théorème de la bijection + Théorème des valeurs intermédiaires|
|:-|
|Soit $f$ une fonction continue et **strictement monotone** sur $I=[a,b]$ de $\mathbb{R}$|
|alors $f$ induit une bijection de $I$ dans $f(I)$.|
|Si de plus $f(a)f(b)<0$|
|alors il existe un **unique** $c \in ]a,b[$ tel que $f(c)=0$.|
|Autrement dit, la fonction $f$ s'annule une seule fois dans $]a,b[$.|

L'étude préalable de la fonction doit donc avoir pour but de séparer les racines en isolant des intervalles sur lesquels la fonction est **strictement monotone** et **change de signe**.

### Méthodes numériques et convergence

Comme évoqué précédemment, l'idée est d'approximer la racine d'une fonction lorsque l'on est incapables de trouver la solution analytiquement.

Les méthodes numériques de recherche de racines sont en général des méthodes **itératives**.

Il s'agit de construire une suite $x_{n+1} = g(x_n)$ telle que $\lim\limits_{n \to \infty} x_n = c$ où $c$ est la racine à approcher.

La performance de ces méthodes est évaluée par leur **vitesse de convergence** :

|Définition|
|:-|
|Soit $c$ la racine recherchée. Posons $e_n = x_n - c$ l'**erreur absolue** à l'itération $n$.|
|La suite est dite **convergente d'ordre $p \geq 1$** si|
|il existe une constante $K>0$ et un indice $n_0 \geq 1$ tels que $\forall n \geq n_0$ :| 
|$\lim\limits_{n \to \infty} \frac{\mid e_{n+1} \mid}{\mid e_n \mid^p} \leq K$|

La convergence est d'autant plus rapide que la valeur de $p$ est grande.
$K$ est le **facteur de convergence** de la suite.

Voici comment on qualifie la convergence en fonction de $$p$$ et $$K$$ :

|Valeur de $p$|Valeur de $K$|Convergence                     |
|:------------|:-----------:|-------------------------------:|
|$1$          |$0$          |Super-linéaire                  |
|$1$          |$]0,1[$      |Linéaire                        |
|$1$          |$1$          |Logarithmique                   |
|$1$          |$>1$         |Sous-linéaire (non-convergence) |
|$2$          |             |Quadratique                     |
|$3$          |             |Cubique                         |
|$4$          |             |Quadratique                     |

Dans la pratique, la racine étant inconnue, **nous ne pouvons pas calculer l'erreur** $e_n$.

C'est pourquoi, à chaque itération, on calcule plutôt le **résidu** $r_n = f(x_n)$.

On considère que la suite est suffisamment proche de la racine si $r_n < \varepsilon$ avec $\varepsilon$ la précision choisie.

La convergence des méthodes itératives de recherche de racines dépend en général du choix de la donnée initiale $x_0$ :

* Une méthode qui converge quelque soit $x_0$ est dite **globalement convergente**.

* Une méthode qui converge seulement lorque $x_0$ est au voisinage de la racine est dite **localement convergente**.

De manière générale, les méthodes localement convergentes on un ordre de convergence plus grand que les méthodes globalement convergentes.

---

## Méthode de la dichotomie

### Algorithme

La méthode de la **dichotomie* est inspirée du théorème des valeurs intermédiaires.
Elle est aussi connue sous le nom de **"méthode de la bissection"**.

Soit $f$ une fonction continue de $[a,b]$ dans $\mathbb{R}$.
On suppose que $f$ admet une unique racine dans $]a,b[$ et que $f(a)f(b)<0$.

Voici l'algorithme sous la forme d'une fonction Python.
Elle prend en entrée :

* `f` la fonction dont on cherche les racines.

* `a` et `b` les bornes de l'intervalle de recherche.

* `n_max` le nombre maximum d'itérations.

* `e` la précision désirée.

On notera `x_n`, `a_n` et `b_n` les variables correspondant à l'estimation de la racine, la borne inférieure et la borne supérieure de l'intervalle de recherche à l'itération `n`.

~~~
def dichotomie(f,a,b,n_max,e):
	#Initialisation de x_n
~~~

### Convergence

### Précision

### Exemple

---

## Méthode de la sécante

---

## Méthode de la fausse position

---

## Méthode de Newton

---

## Méthode du point fixe