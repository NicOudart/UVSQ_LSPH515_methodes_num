import numpy as np

pos_rec = np.array([4205,158,4777],dtype=np.float64) #Coordonnées du récepteur GPS

pos_sat1 = np.array([14000,4000,25000],dtype=np.float64) #Coordonnées du satellite GPS 1
pos_sat2 = np.array([24000,6000,15000],dtype=np.float64) #Coordonnées du satellite GPS 2
pos_sat3 = np.array([9000,-14000,21000],dtype=np.float64) #Coordonnées du satellite GPS 3
pos_sat4 = np.array([10000,16000,19000],dtype=np.float64) #Coordonnées du satellite GPS 4

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

#Calcul du déterminant---------------------------------------------------------

def det_3(A):
    
    return A[0,0]*A[1,1]*A[2,2]+A[0,1]*A[1,2]*A[2,0]+A[0,2]*A[1,0]*A[2,1]-A[0,2]*A[1,1]*A[2,0]-A[0,1]*A[1,0]*A[2,2]-A[0,0]*A[1,2]*A[2,1]

#Méthode de Cramer-------------------------------------------------------------

def cramer_3(A,b):
    
    #Vérification des dimensions de A (3x3) et b (3) :
    if (np.shape(A)!=(3,3))or(len(b)!=3):
        
        raise ValueError("Le système n'est pas de Cramer")
    
    #Calculer le déterminant de A :
    det_A = det_3(A)
    
    #Vérifier que le système admet bien une unique solution :
    if det_A==0:
        
        raise ValueError("Le système n'admet pas une solution unique !")
        
    #Initialiser le vecteur qui contiendra la solution :
    x = np.array([0,0,0],dtype=np.float64)
    
    #Boucle sur les 3 colonnes de la matrice A :
    for i in range(3):
        
        #Remplir la matrice A_i avec les éléments de A :
        A_i = np.copy(A)
        
        #Remplacer la i-ème colonne de A_i avec les éléments de b :
        A_i[:,i] = b
        
        #
        x[i] = det_3(A_i)/det_A
    
    return x

print('Cramer method solution = '+str(cramer_3(A,b)))