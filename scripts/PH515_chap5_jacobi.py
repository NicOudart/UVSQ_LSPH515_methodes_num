import numpy as np

pos_rec = np.array([4205,158,4777],dtype=np.float64) #Coordonnées du récepteur GPS

pos_sat1 = np.array([15000,13000,18000],dtype=np.float64) #Coordonnées du satellite GPS 1
pos_sat2 = np.array([1000,6000,24000],dtype=np.float64) #Coordonnées du satellite GPS 2
pos_sat3 = np.array([19000,2000,19000],dtype=np.float64) #Coordonnées du satellite GPS 3
pos_sat4 = np.array([12000,12000,26000],dtype=np.float64) #Coordonnées du satellite GPS 4

#Génération des temps de retard------------------------------------------------

t1 = ((sum((pos_rec-pos_sat1)**2))**0.5)/3e5 #Temps de retard associé au satellite GPS 1
t2 = ((sum((pos_rec-pos_sat2)**2))**0.5)/3e5 #Temps de retard associé au satellite GPS 2
t3 = ((sum((pos_rec-pos_sat3)**2))**0.5)/3e5 #Temps de retard associé au satellite GPS 3
t4 = ((sum((pos_rec-pos_sat4)**2))**0.5)/3e5 #Temps de retard associé au satellite GPS 4


#Génération des matrices A et b------------------------------------------------

A = np.vstack((pos_sat2-pos_sat1,pos_sat3-pos_sat1,pos_sat4-pos_sat1))

b_row1 = 0.5*(sum(pos_sat2**2)-sum(pos_sat1**2)-(3e5*t2)**2+(3e5*t1)**2)
b_row2 = 0.5*(sum(pos_sat3**2)-sum(pos_sat1**2)-(3e5*t3)**2+(3e5*t1)**2)
b_row3 = 0.5*(sum(pos_sat4**2)-sum(pos_sat1**2)-(3e5*t4)**2+(3e5*t1)**2)

b = np.array([b_row1,b_row2,b_row3])

#Méthode de Jacobi-------------------------------------------------------------

def jacobi(A,b,x_0,n_max,e):
    
    #Récupérer les dimensions de la matrice A :
    m,n = np.shape(A)
    
    #Vérification des dimensions de A (nxn) et b (n) :
    if (m!=n)or(len(b)!=n):
        
        raise ValueError("Le système n'est pas de Cramer")
    
    #Récupérer les valeurs sur la diagonale de A :
    A_diag = np.diag(A)
    
    #Créer la matrice diagonale D, ayant les mêmes valeurs que la diagonale de A :
    D = np.diag(A_diag)
    
    #Créer la matrice E+F, ayant des zéros sur sa diagonale et les mêmes valeurs
    #que A partout ailleurs :
    EF = A-D
    
    #Vérifier si la matrice A est à diagonale strictement dominante :
    for i in range(n):
        
        if sum(abs(EF[:,i]))>=abs(A_diag[i]):
            
            print("Attention : la matrice A n'est pas à diagonale strictement dominante")
            break
        
    #Initialisation des variables :
    n = 0 #Nombre d'itérations
    x_n_old = np.copy(x_0) #Estimation de la solution à l'itération n-1
    x_n = (b-np.dot(EF,x_n_old))/A_diag #Estimation de la solution à l'itération n
    r_n = np.dot(A,x_n)-b #Résidu
        
    #Itérations de l'algorithme de Jacobi
	#tant qu'une des conditions d'arrêt n'est pas atteinte :
    while (n<n_max)and(np.linalg.norm(x_n-x_n_old,ord=2)>e)and(np.linalg.norm(r_n,ord=2)>e):
        
        #Mettre à jour l'estimation de la solution :
        x_n_old = np.copy(x_n) #Itération n
        x_n = (b-np.dot(EF,x_n_old))/A_diag #Iteration n+1
        
        #Incrémenter le nombre d'itérations :
        n+=1
        
        #Mettre à jour le résidu :
        r_n = np.dot(A,x_n)-b
                    
    #Renvoyer l'estimation de la solution du système et le résidu :
    return x_n,r_n

#Application-------------------------------------------------------------------

x_0 = np.array([0,0,0],dtype=np.float64)

x_n,r_n = jacobi(A,b,x_0,100,1e-3)