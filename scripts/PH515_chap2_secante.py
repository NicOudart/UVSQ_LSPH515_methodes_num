def f(x):
    
    return (x**2)-2

def secante(f,a,b,n_max,e):

	#Initialisation des variables :
	n = 1 #Nombre d'itérations
	x_n = a #Estimation de la racin à l'itération n
	x_n_old = b #Estimation de la racin à l'itération n-1
	r_n = f(x_n) #Résidu
	
	#Itérations de l'algorithme de la sécante
	#tant qu'une des conditions d'arrêt n'est pas atteinte :
	while (n<n_max)and(abs(x_n-x_n_old)>e)and(abs(r_n)>e):
	
		#Calculer la pente : 
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

x_n,r_n = secante(f,0,2,100,1e-6)