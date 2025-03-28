def f(x):
    
    return (x**2)-2

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

x_n,r_n = fausse_position(f,0,2,100,1e-6)