# Chapitre I : Introduction aux méthodes numériques

Ce chapitre est une introduction aux enjeux des méthodes numériques.

---

Nous allons détailler dans ce chapitre  les différentes étapes de résolution d'un problème en Physique, avec un exemple simple.

## Le problème

_Un groupe d'étudiants de l'UVSQ réalise un projet de ballon sonde pour mesurer la vitesse du vent dans la stratosphère au cours du temps_.
_On considère qu'après une phase d'ascension, leur ballon a atteint une altitude stable dans la stratosphère, et que son déplacement est uniquement lié au vent_.
_La nacelle de leur ballon contient une balise GPS, qui leur permet de mesurer le déplacement du ballon au cours du temps_.

Ceci est le **problème physique** auquel les étudiants sont confrontés.

## Modélisation du problème

La première étape est de **traduire ce problème physique en un modèle** par le biais d'équations mathématiques : 

_Les étudiants font l'hypothèse que leur ballon stratosphérique se déplace à la vitesse du vent_.

Notons $p(t)$ la fonction associant à chaque instant $t$ en [s] à partir du début de la mesure le déplacement du ballon en [m].

La traduction mathématique du problème est donc que la vitesse du vent au cours du temps $W(t)$ est liée à $p(t)$ par l'équation :

$W(t) = \frac{d}{dt} p(t)$

Comme la fonction $p(t)$ n'a pas d'expression analytique connue (il s'agit d'une mesure d'un capteur), les étudiants devront utiliser une méthode numérique pour estimer la dérivée de $p(t)$, et donc la solution du problème.

## Discrétisation du problème

Un ordinateur ne pouvant gérer des objets continus, l'application d'une méthode numérique nécessite une **discrétisation**.

Dans le cas des étudiants de l'UVSQ, la discrétisation est réalisée au moment de l'échantillonnage :

_Les étudiants décident d'enregistrer une mesure de déplacement toutes les 10 s._

Ils n'ont donc pas accès à toutes les valeurs possibles de $p(t)$, mais à des valeurs discrètes régulièrement espacées $p(t_i)$ avec $t_i = 0, 10, 20, 30, ... 10 \times (N-1)$ et $N$ le nombre d'échantillons.

De même, ils n'estimeront pas toutes les valeurs de $W(t)$, mais des valeurs discrètes $W(t_i)$.

On notera $h = 10 s$ le pas de discrétisation du problème.

## Choix d'une méthode numérique

Nos étudiants de l'UVSQ on besoin ici d'une **méthode numérique** d'estimation de la dérivée d'une fonction.

Pour approcher la dérivée  décident d'employer une méthode de "_différences décentrées à droite_" :

$W(t_i) \approx \frac{p(t_{i+1})-p(t_i)}{h}$ avec $t_{i+1} = t_i + h$ et $i = 0, 1, ..., N-2$

Lorsque l'on choisi une méthode numérique pour répondre à un problème, il convient de se poser les questions suivantes :

### La méthode est-elle **applicable** ? 

Dans notre exemple, il faut que la fonction $p(t)$ soit doublement dérivable sur chaque intervalle $[t_i,t_{i+1}]$.

### La solution est-elle **unique** ? 

Dans notre exemple, l'unicité est évidente, mais ce n'est pas toujours le cas.

### Quelles sont les sources d'erreur, et quelle est l'erreur attendue ?

### La méthode **converge**-t-elle vers la solution ? Avec quelle vitesse ?

Une méthode numérique est dites **convergente** si l'écart entre la solution approchée et la solution exacte tend vers 0 quand le pas de discrétisation $h$ tend vers 0.

Si de plus, l'erreur absolue 

Dans notre exemple, on peut montrer que pour chaque intervalle $[t_i,t_{i+1}]$, l'erreur est majorée par $\frac{h}{2} sup_{t \in [t_i,t_{i+1}]} \mid \frac{d^2}{dt^2} p(t)\mid$. 
La méthode converge donc vers la solution lorsque que le pas de discrétisation $h$ diminue.
Nous verrons dans la suite que cette convergence est dite "d'ordre 1".

### La méthode est-elle stable ?

### Quelles sont les demandes de ressources informatique de la méthode : temps de calcul et mémoire ?

## Implémentation de l'algorithme

## Execution et analyse du résultat