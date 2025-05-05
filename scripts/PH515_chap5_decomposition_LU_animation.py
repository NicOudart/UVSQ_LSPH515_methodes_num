import numpy as np

pos_rec = np.array([4205,158,4777],dtype=np.float64) #Coordonnées du récepteur GPS

pos_sat1 = np.array([14000,4000,25000],dtype=np.float64) #Coordonnées du satellite GPS 1
pos_sat2 = np.array([9000,-14000,21000],dtype=np.float64) #Coordonnées du satellite GPS 2
pos_sat3 = np.array([24000,6000,15000],dtype=np.float64) #Coordonnées du satellite GPS 3
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

#Decomposition LU--------------------------------------------------------------
    
#Récupérer les dimensions de la matrice A :
m,n = np.shape(A)

#Vérification des dimensions de A (nxn) :
if (m!=n):
    
    raise ValueError("La matrice A n'est pas carrée")
 
#Copier A pour ne pas modifier la matrice originale :
U = np.copy(A)

#Initialiser la matrice L comme la matrice identité (nxn) :
L = np.eye(n)

print('Système initial :')
print('L = '+str(L))
print('U = '+str(U))
print('-----------')

#Boucle sur les colonnes de la matrice A, jusqu'à l'avant-dernière :
for j in range(n-1):
    
    print('Itération colonne '+str(j+1))
    
    #Sélection du pivot comme étant la valeur sur la diagonale de la j-ème colonne :
    pivot = U[j,j]
    print('Pivot = '+str(pivot))
    
    #On vérifie que le pivot n'est pas nul :
    if pivot!=0:
        
        #Boucle sur les lignes sous le pivot :
        for k in range(j+1,n):
            
            #Sauvegarde du coefficient d'élimination de Gauss dans L :
            L[k,j] = U[k,j]/pivot
            print('Ajouter '+str(U[k,j])+'/'+str(pivot)+' à la ligne '+str(k+1)+' de L')
            
            #Opérations d'élimination de Gauss sur les lignes de A en 
            #utilisant le pivot :
            print('Opération L'+str(k+1)+' = L'+str(k+1)+' - L'+str(j+1)+'x'+str(U[k,j])+'/'+str(pivot)+' sur U')
            U[k,:] = U[k,:] - U[j,:]*L[k,j]
            
    print('L = '+str(L))
    print('U = '+str(U))
    print('-----------')