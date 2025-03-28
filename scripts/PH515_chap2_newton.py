def f(x):

	return (x**2)-2

def f_derivee(x):

	return 2*x

def newton(f,f_derivee,x_0,n_max,e):

	#Vérifier que le point de départ de la recherche est possible :
    if f_derivee(x_0)==0:
        raise ValueError("Mauvaise initialisation, f'(x_0) = 0")

	#Initialisation des variables :
    n = 0 #Nombre d'itérations
    x_n_old = x_0 #Estimation de la racine à l'itération n-1
    x_n = x_n_old-f(x_n_old)/f_derivee(x_n_old) #Estimation de la racine à l'itération n
    r_n = f(x_n) #Résidu
	
	#Itérations de l'algorithme de Newton
	#tant qu'une des conditions d'arrêt n'est pas atteinte :
    while (n<n_max)and(abs(x_n-x_n_old)>e)and(abs(r_n)>e):
	
		#Mettre à jour l'estimation de la racine :
        x_n_old = x_n #Itération n
        x_n = x_n-f(x_n)/f_derivee(x_n) #Iteration n+1
		
		#Incrémenter le nombre d'itérations :
        n+=1
		
		#Mettre à jour le résidu :
        r_n = f(x_n)

	#Renvoyer l'estimation de la racine et le résidu :
    return x_n,r_n

x_n,r_n = newton(f,f_derivee,2,100,1e-6)