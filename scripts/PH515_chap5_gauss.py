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

#Gauss sans pivot--------------------------------------------------------------

def gauss_sans_pivot(A,b):
    
    #Récupérer les dimensions de la matrice A :
    m,n = np.shape(A)
    
    #Vérification des dimensions de A (nxn) et b (n) :
    if (m!=n)or(len(b)!=n):
        
        raise ValueError("Le système n'est pas de Cramer")
     
    #Copier A et b pour ne pas modifier les matrices originales :
    A_2 = np.copy(A)
    b_2 = np.copy(b)
    
    #Boucle sur les colonnes de la matrice A, jusqu'à l'avant-dernière :
    for j in range(n-1):
        
        #Sélection du pivot comme étant la valeur sur la diagonale de la j-ème colonne :
        pivot = A_2[j,j]
        
        #On vérifie que le pivot n'est pas nul :
        if pivot!=0:
            
            #Boucle sur les lignes sous le pivot :
            for k in range(j+1,n):
                
                #Opérations sur les lignes de A et b en utilisant le pivot :
                b_2[k] = b_2[k] - b_2[j]*A_2[k,j]/pivot
                A_2[k,:] = A_2[k,:] - A_2[j,:]*A_2[k,j]/pivot
    
    #Renvoyer les matrices A et b modifiées :
    return A_2,b_2

#Gauss pivot partiel-----------------------------------------------------------

def gauss_pivot_partiel(A,b):
    
    #Récupérer les dimensions de la matrice A :
    m,n = np.shape(A)
    
    #Vérification des dimensions de A (nxn) et b (n) :
    if (m!=n)or(len(b)!=n):
        
        raise ValueError("Le système n'est pas de Cramer")
     
    #Copier A et b pour ne pas modifier les matrices originales :
    A_2 = np.copy(A)
    b_2 = np.copy(b)
    
    #Boucle sur les colonnes de la matrice A, jusqu'à l'avant-dernière :
    for j in range(n-1):
        
        #Sélection du pivot comme étant la valeur maximale en absolu sur la colonne, 
        #sur la j-ième ligne ou en dessous :
        idx_pivot = np.argmax(abs(A_2[j:,j]))+j #Indice de la ligne du pivot
        pivot = A_2[idx_pivot,j] #Valeur du pivot
        
        #Si le pivot n'est pas sur la j-ième ligne, échanger la j-ième et la
        #ligne du pivot :
        if idx_pivot!=j:
            A_2[[j,idx_pivot]] = A_2[[idx_pivot,j]] #Pour la matrice A
            b_2[[j,idx_pivot]] = b_2[[idx_pivot,j]] #Pour le vecteur b
        
        #On vérifie que le pivot n'est pas nul :
        if pivot!=0:
            
            #Boucle sur les lignes sous le pivot :
            for k in range(j+1,n):
                
                #Opérations sur les lignes de A et b en utilisant le pivot :
                b_2[k] = b_2[k] - b_2[j]*A_2[k,j]/pivot
                A_2[k,:] = A_2[k,:] - A_2[j,:]*A_2[k,j]/pivot
    
    #Renvoyer les matrices A et b modifiées :
    return A_2,b_2

#Gauss pivot total-------------------------------------------------------------

def gauss_pivot_total(A,b):
    
    #Récupérer les dimensions de la matrice A :
    m,n = np.shape(A)
    
    #Vérification des dimensions de A (nxn) et b (n) :
    if (m!=n)or(len(b)!=n):
        
        raise ValueError("Le système n'est pas de Cramer")
     
    #Copier A et b pour ne pas modifier les matrices originales :
    A_2 = np.copy(A)
    b_2 = np.copy(b)
    
    idx_x = np.arange(n)
    
    #Boucle sur les colonnes de la matrice A, jusqu'à l'avant-dernière :
    for j in range(n-1):
        
        #Sélection du pivot comme étant la valeur maximale en absolu sur la 
        #portion de matrice non-triangularisée :
        idx_pivot = np.argmax(abs(A_2[j:,j:])) #Indice du pivot
        ligne_pivot = idx_pivot//(n-j)+j
        colonne_pivot = idx_pivot%(n-j)+j
        pivot = A_2[ligne_pivot,colonne_pivot] #Valeur du pivot
        
        #Si le pivot n'est pas sur la j-ième ligne, échanger la j-ième et la
        #ligne du pivot :
        if ligne_pivot!=j:
            A_2[[j,ligne_pivot]] = A_2[[ligne_pivot,j]] #Pour la matrice A
            b_2[[j,ligne_pivot]] = b_2[[ligne_pivot,j]] #Pour le vecteur b
            
        #Si le pivot n'est pas sur la j-ième colonne, échanger la j-ième et la
        #colonne du pivot :
        if colonne_pivot!=j:
            A_2[:,[j,colonne_pivot]] = A_2[:,[colonne_pivot,j]] #Pour la matrice A
            idx_x[[j,colonne_pivot]] = idx_x[[colonne_pivot,j]] #Pour les éléments de x
        
        #On vérifie que le pivot n'est pas nul :
        if pivot!=0:
            
            #Boucle sur les lignes sous le pivot :
            for k in range(j+1,n):
                
                #Opérations sur les lignes de A et b en utilisant le pivot :
                b_2[k] = b_2[k] - b_2[j]*A_2[k,j]/pivot
                A_2[k,:] = A_2[k,:] - A_2[j,:]*A_2[k,j]/pivot
    
    #Renvoyer les matrices A et b modifiées :
    return A_2,b_2,idx_x

#Algorithme de remontée--------------------------------------------------------

def remontee(A,b):
    
    #Récupérer le nombre n d'équations / inconnues du système :
    n = len(A)
    
    #Initialiser le vecteur qui contiendra les solutions du système :
    x = np.zeros(n,dtype=np.float64)
    
    #Boucle sur les lignes de la matrice A, de n à 1 :
    for i in range(n-1,-1,-1):
        
        #Détermination de la i-ème inconnue :
        x[i] = (b[i]-sum(A[i,i+1:n]*x[i+1:n]))/A[i,i]
    
    #Renvoyer le vecteur contenant les solutions du système :
    return x

#Applications------------------------------------------------------------------

A_1,b_1 = gauss_sans_pivot(A,b)
x_1 = remontee(A_1,b_1)

print('Gauss sans pivot :')
print('A = '+str(A_1))
print('b = '+str(b_1))
print('x = '+str(x_1))

A_2,b_2 = gauss_pivot_partiel(A,b)
x_2 = remontee(A_2,b_2)

print('Gauss avec pivot partiel :')
print('A = '+str(A_2))
print('b = '+str(b_2))
print('x = '+str(x_2))

A_3,b_3,idx_x = gauss_pivot_total(A,b)
x_3 = remontee(A_3,b_3)

print('Gauss avec pivot total :')
print('A = '+str(A_3))
print('b = '+str(b_3))
print('x = '+str(x_3))