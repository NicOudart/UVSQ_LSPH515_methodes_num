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

**Ce chapitre présentera un panel de méthodes numériques de recherche de racines**, que nous appliquerons à un même exemple pour illustration.

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

Voici comment on qualifie la convergence en fonction de $p$ et $K$ :

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

### Exemple de problème

Au cours de ce chapitre, nous appliquerons les différentes méthodes numériques de recherche de racines à un même exemple : **l'estimation de $\sqrt{2}$**.

$\sqrt{2}$ est un nombre irrationel, dont l'approximation est un problème depuis l'antiquité, notamment parce qu'il correspond à l'hypothénuse d'un carré de côté 1.

Par définition, $\sqrt{2}$ et $-\sqrt{2}$ sont les solutions de l'équation $x^2 = 2$.

Résoudre cette équation revient à résoudre $x^2 - 2 = 0$.

Pour calculer une approximation de $\sqrt{2}$, on peut donc chercher les racines de la fonction **$f(x) = x^2 - 2$** de $\mathbb{R}$ dans $\mathbb{R}$.

![Graphique de f](img/Chap2_exemple_fonction.png)

La voici sous la forme d'une fonction Python :

~~~
def f(x):

	return (x**2)-2
~~~

Cette fonction est continue et dérivable sur $\mathbb{R}$, et sa dérivée est $f'(x) = 2x$. On peut en déduire ses variations :

![Tableau de variations de f](img/Chap2_exemple_tab_var.png)

* Théorème de la bijection : $f$ est continue et strictement monotone sur $I=[1,2]$, donc $f$ induit une bijection de $I$ dans $f(I)$.

* Théorème des valeurs intermédiaires : $f(1)=-1$ et $f(2)=2$, d'où $f(1)f(2)<0$.

Donc, il existe une seule racine de $f$ dans $]1,2[$, et nous savons que cette racine est $\sqrt{2}$.

C'est pourquoi dans la suite de ce chapitre, sauf indication contraire, nous chercherons la racine de $f$ se trouvant sur l'intervalle $]1,2[$.

---

## Méthode de la dichotomie

### Algorithme

La méthode de la **dichotomie** est inspirée du théorème des valeurs intermédiaires.
Elle est aussi connue sous le nom de **"méthode de la bissection"**.

Soit $f$ une fonction continue de $[a,b]$ dans $\mathbb{R}$.
On suppose que $f$ admet une unique racine dans $]a,b[$ et que $f(a)f(b)<0$.

Voici l'algorithme sous la forme d'une fonction Python.

Elle prend en entrée :

* `f` la fonction dont on cherche les racines.

* `a` et `b` les bornes de l'intervalle de recherche.

* `n_max` le nombre maximum d'itérations.

* `e` la précision désirée.

On notera les variables à l'itération `n` : 

* `x_n` l'estimation de la racine.

* `a_n` et `b_n` les bornes de l'intervalle de recherche.

* `r_n` le résidu.

~~~
def dichotomie(f,a,b,n_max,e):

	#Initialisation des variables :
	n = 0 #Nombre d'itérations
	x_n = (a+b)/2 #Estimation de la racine
	a_n = a #Borne inférieure de l'intervalle de recherche
	b_n = b #Borne supérieure de l'intervalle de recherche
	r_n = f(x_n) #Résidu
	
	#Itérations de l'algorithme de dichotomie
	#tant qu'une des conditions d'arrêt n'est pas atteinte :
	while (n<n_max)and(abs(a_n-b_n)>e)and(abs(r_n)>e):
	
		#Si la racine est dans ]a_n,x_n[ alors on remplace b_n par x_n:
		if (f(a_n)*f(x_n))<0:
			b_n = x_n
		
		#Si la racine est dans ]x_n,b_n[ alors on remplace a_n par x_n:
		if (f(x_n)*f(b_n))<0:
			a_n = x_n
		
		#Incrémenter le nombre d'itérations :
		n+=1
		
		#Mettre à jour l'estimation de la racine et le résidu :
		x_n = (a_n+b_n)/2 
		r_n = f(x_n)
	
	#Renvoyer l'estimation de la racine et le résidu :
	return x_n,r_n
~~~

### Convergence

La longueur de l'intervalle de recherche est divisée par 2 à chaque itération :

$I_n = \mid b_n-a_n \mid = \frac{\mid b-a \mid}{2^n}$

Donc, l'erreur absolue à l'itération $n \geq 0$, $e_n=x_n-c$

$\mid e_n \mid \leq \frac{I_n}{2} = \frac{\mid b-a \mid}{2^{n+1}}$

Ce qui entraine : $\lim\limits_{n \to \infty} \mid e_n \mid = 0$
et $\frac{\mid e_{n+1} \mid}{\mid e_n \mid} \leq \frac{1}{2}$

La méthode de la dichotomie est donc **globalement convergente** : elle converge quelque soit le point de départ.
(Si la $f$ a plusieurs racines dans $[a,b]$, la racine trouvée dépendra de l'intervalle).

Sa convergence est **linéaire**, donc elle est relativement **lente**.
C'est pourquoi on utilise souvent cette méthode juste pour initialiser une méthode plus rapide.

### Précision

En fonction de la précision souhaitée $\varepsilon$, on peut calculer le nombre d'itérations $m$ pour approcher la racine.

On cherche $m \in \mathbb{N}$ tel que :
$\mid e_{m-1} \mid \leq \frac{\mid b-a \mid}{2^m} \leq \varepsilon$

Donc tel que :
$2^m \geq \frac{\mid b-a \mid}{\varepsilon}$

Soit : 
$m \geq log_2{\frac{\mid b-a \mid}{\varepsilon}} = \frac{ln{\frac{\mid b-a \mid}{\varepsilon}}}{ln(2)} \approx 1.4427 ln{\frac{\mid b-a \mid}{\varepsilon}}$

### Exemple

Voici les 4 premières itérations de la méthode de la dichotomie appliquée à notre problème exemple, pour un intervalle initial $[1,2]$ :

![Exemple d'application de la dichotomie](img/Chap2_exemple_dichotomie.gif)

**Exercice :**

En adaptant la fonction Python donnée précédemment pour la méthode de la dichotomie, avec un intervalle initial $[1,2]$, estimez la valeur de $\sqrt{2}$ avec une précision de $10^{-6}$.
Combien d'itérations sont nécessaires pour obtenir cette précision ? Retrouvez-vous bien le nombre d'itérations théorique ?

---

## Avant-propos : les méthodes linéarisées

Les méthodes qui seront présentées dans la suite de ce chapitre sont des **méthodes linéarisées**.
Ce type de méthodes s'appuit sur le dévloppement de Taylor de $f$ autour de sa racine $c$ :

$f(c) = 0 = f(x') + (c-x')f'(\xi)$ où $\xi \in [x',c]$ et $f(x') \neq 0$

D'où $c = x' - \frac{f(x')}{f'(\xi)}$

Donc, si on connait $\xi$, on peut déterminer $c$ à partir de $x'$.

D'un point de vue géométrique, la racine $c$ est à l'intersection entre la droite passant par le point $(x',f'(x'))$ et de pente $f'(\xi)$ et donc d'équation :

$y = f'(\xi) x + f(x') - f'(\xi) x'$

Et l'axe $(Ox)$ donc d'équation $y = 0$.

![Illustration des méthodes linéarisées](img/Chap2_methodes_linearisees.png){width="350"}

D'où la méthode itérative suivante :

$f(x_n) + (x_{n+1} - x_n) q_n = 0$

ou encore l'**équation de récurrence** :

$x_{n+1} = x_n - \frac{f(x_n)}{q_n}$

où $q_n$ est une approximation de $f'(\xi)$.

L'idée des méthodes linéarisées est donc :

* D'approcher une fonction non-linéaire par une droite.

* Déterminer à chaque itération $x_{n+1}$ comme l'intersection entre l'axe $(Ox)$ et la droite de pente $q_n$ passant par le point $(x_n,f(x_n))$.

Les méthodes linéarisées (méthode de la sécante, méthode de la fausse position, méthode de Newton, etc.) se différentient par **le choix de $q_n$**.

---

## Méthode de la sécante

### Algorithme

La **méthode de la sécante** est une méthode linéarisée pour laquelle :

$q_n = \frac{f(x_n)-f(x_{n-1})}{x_n-x_{n-1}}$

Cette suite correspond à la droite passant par les points $(x_n,f(x_n))$ et $(x_{n-1},f(x_{n-1}))$.

Soit $f$ une fonction continue de $[a,b]$ dans $\mathbb{R}$.
On suppose que $f$ admet une unique racine dans $]a,b[$ et que $f(a)f(b)<0$.

Voici l'algorithme sous la forme d'une fonction Python.

Elle prend en entrée :

* `f` la fonction dont on cherche les racines.

* `a` et `b` les bornes de l'intervalle de recherche.

* `n_max` le nombre maximum d'itérations.

* `e` la précision désirée.

On notera les variables à l'itération `n` : 

* `x_n` l'estimation de la racine.

* `a_n` et `b_n` les bornes de l'intervalle de recherche.

* `r_n` le résidu.

~~~
def secante(f,a,b,n_max,e):

	#Initialisation des variables :
	n = 1 #Nombre d'itérations
	x_n = a #Estimation de la racin à l'itération n
	x_n_old = b #Estimation de la racin à l'itération n-1
	r_n = f(x_n) #Résidu
	
	#Itérations de l'algorithme de la sécante
	#tant qu'une des conditions d'arrêt n'est pas atteinte :
	while (n<n_max)and(abs(x_n-x_n_old)>e)and(abs(r_n)>e):
	
		#Calculer la pente de la droite : 
		q_n = (f(x_n)-f(x_n_old))/(x_n-x_n_old)
		
		#Mettre à jour l'estimation de la racine :
		x_n_old = x_n #Iteration n
		x_n = x_n - f(x_n)/q_n #Iteration n+1
		
		#Incrémenter le nombre d'itérations :
		n+=1
		
		#Mettre à jour le résidu :
		r_n = f(x_n)
		
	#Renvoyer l'estimation de la racine et le résidu :
	return x_n,r_n
~~~

### Convergence

La méthode de la sécante est **convergente localement**.

Si les données initiales sont assez proches de la racine $c$, et que $f'(c) \neq 0$, alors on peut démontrer qu'elle converge avec un ordre $p = \frac{1+\sqrt(5)}{2}$.
Cette valeur est connue sous le nom de "nombre d'or".

### Exemple

Voici les 5 premières itérations de la méthode de la sécante appliquée à notre problème exemple.
L'intervalle initial est ici de [0,2] pour des raisons de lisibilité :

![Exemple d'application de la sécante](img/Chap2_exemple_secante.gif)

**Exercice :**

En adaptant la fonction Python donnée précédemment pour la méthode de la sécante, avec un intervalle initial $[1,2]$, estimez la valeur de $\sqrt{2}$ avec une précision de $10^{-6}$.
Combien d'itérations sont nécessaires pour obtenir cette précision ? Comparez cette valeur à celle obtenue pour la méthode de la dichotomie.

---

## Méthode de la fausse position

### Algorithme

La **méthode de la fausse position** est une méthode linéarisée pour laquelle :

$q_n = \frac{f(x_n)-f(x_n')}{x_n-x_n'}$

Cette suite correspond à la droite passant par les points $(x_n,f(x_n))$ et $(x_n',f(x_n'))$, où $n'$ est le plus grand indice inférieur à $n$ tel que $f(x_n)f(x_n')<0$.

Il s'agit d'un mélange entre la méthode de la dichotomie et la méthode de la sécante.
On l'appelle aussi **méthode de regula falsi** ou **méthode de Lagrange**.

Soit $f$ une fonction continue de $[a,b]$ dans $\mathbb{R}$.
On suppose que $f$ admet une unique racine dans $]a,b[$ et que $f(a)f(b)<0$.

Voici l'algorithme sous la forme d'une fonction Python.

Elle prend en entrée :

* `f` la fonction dont on cherche les racines.

* `a` et `b` les bornes de l'intervalle de recherche.

* `n_max` le nombre maximum d'itérations.

* `e` la précision désirée.

On notera les variables à l'itération `n` : 

* `x_n` l'estimation de la racine.

* `a_n` et `b_n` les bornes de l'intervalle de recherche.

* `r_n` le résidu.

~~~
def fausse_position(f,a,b,n_max,e):

	#Initialisation des variables :
	n = 0 #Nombre d'itérations
	a_n = a #Borne inférieure de l'intervalle de recherche
	b_n = b #Borne supérieure de l'intervalle de recherche
	q_n = (f(b_n)-f(a_n))/(b_n-a_n) #Initialiser la pente de la droite
	x_n = a_n-f(a_n)/q_n #Estimation de la racine
	r_n = f(x_n) #Résidu
	
	#Itérations de l'algorithme de la fausse-position
	#tant qu'une des conditions d'arrêt n'est pas atteinte :
	while (n<n_max)and(abs(a_n-b_n)>e)and(abs(r_n)>e):
	
		#Si la racine est dans ]a_n,x_n[ alors on remplace b_n par x_n:
		if (f(a_n)*f(x_n))<0:
			b_n = x_n
		
		#Si la racine est dans ]x_n,b_n[ alors on remplace a_n par x_n:
		if (f(x_n)*f(b_n))<0:
			a_n = x_n
	
		#Incrémenter le nombre d'itérations :
		n+=1
		
		#Calcul de la nouvelle pente de la droite :
		q_n = (f(b_n)-f(a_n))/(b_n-a_n)
	
		#Mettre à jour l'estimation de la racine et le résidu :
		x_n = x_n-f(x_n)/q_n
		r_n = f(x_n)
	
	#Renvoyer l'estimation de la racine et le résidu :
	return x_n,r_n
~~~

---

## Méthode de Newton

---

## Méthode du point fixe