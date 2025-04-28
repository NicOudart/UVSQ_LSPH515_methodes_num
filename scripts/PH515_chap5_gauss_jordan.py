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

#Gauss-Jordan pivot partiel----------------------------------------------------

def gauss_jordan(A,b):
    
    #Récupérer les dimensions de la matrice A :
    m,n = np.shape(A)
    
    #Vérification des dimensions de A (nxn) et b (n) :
    if (m!=n)or(len(b)!=n):
        
        raise ValueError("Le système n'est pas de Cramer")
     
    #Copier A et b pour ne pas modifier les matrices originales :
    A_2 = np.copy(A)
    b_2 = np.copy(b)
    
    #Boucle sur les colonnes de la matrice A :
    for j in range(n):
        
        #Sélection du pivot comme étant la valeur maximale en absolu sur la colonne, 
        #sur la j-ième ligne ou en dessous :
        idx_pivot = np.argmax(abs(A_2[j:,j]))+j #Indice de la ligne du pivot
        pivot = A_2[idx_pivot,j] #Valeur du pivot
        
        #On vérifie que le pivot n'est pas nul :
        if pivot!=0:
        
            #Division de la ligne du pivot par le pivot, pour que le pivot soit égal à 1 :
            A_2[idx_pivot,:] = A_2[idx_pivot,:]/pivot #Pour la matrice A
            b_2[idx_pivot] = b_2[idx_pivot]/pivot #Pour le vecteur b
            
            #Si le pivot n'est pas sur la j-ième ligne, échanger la j-ième et la
            #ligne du pivot :
            if idx_pivot!=j:
                A_2[[j,idx_pivot]] = A_2[[idx_pivot,j]] #Pour la matrice A
                b_2[[j,idx_pivot]] = b_2[[idx_pivot,j]] #Pour le vecteur b
                
            #Boucle sur toutes les lignes sauf celle du pivot :
            for k in range(n):
                if k!=j:
                
                    #Opérations sur les lignes de A et b :
                    b_2[k] = b_2[k] - b_2[j]*A_2[k,j]
                    A_2[k,:] = A_2[k,:] - A_2[j,:]*A_2[k,j]
    
    #Renvoyer les matrices A et b modifiées :
    return A_2,b_2

#Applications------------------------------------------------------------------

A_2,b_2 = gauss_jordan(A,b)

print('Gauss avec pivot partiel :')
print('A = '+str(A_2))
print('b = '+str(b_2))