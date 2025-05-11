import numpy as np
import matplotlib.pyplot as plt

cas = 2 #Positions possibles du satellite GPS 1 présentées dabs le cours

if cas==0:
    pos_sat1 = np.array([14000,4000,25000],dtype=np.float64) #Coordonnées du satellite GPS 1

if cas==1:
    pos_sat1 = np.array([15000,13000,18000],dtype=np.float64) #Coordonnées du satellite GPS 1

if cas==2:
    pos_sat1 = np.array([23000,7000,20000],dtype=np.float64) #Coordonnées du satellite GPS 1

pos_rec = np.array([4205,158,4777],dtype=np.float64) #Coordonnées du récepteur GPS

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

#Choix du paramètre de relaxation----------------------------------------------

D = np.diag(np.diag(A))
E = D-np.tril(A)
F = D-np.triu(A)

omega = np.linspace(0.001,2,2000)
rho = np.zeros(2000)

for i in range(2000):
    
    C = np.dot(np.linalg.inv(D-omega[i]*E),(1-omega[i])*D+omega[i]*F)
    
    rho[i] = np.max(abs(np.linalg.eigvals(C)))
    
idx_min = np.argmin(rho)

plt.fill([0,0,2,2],[1,np.max(rho),np.max(rho),1],'pink')    
plt.plot(omega,rho,'r-')
plt.plot([omega[idx_min],omega[idx_min]],[0,rho[idx_min]],'r--')
plt.scatter(omega[idx_min],rho[idx_min],color='r',marker='o')
plt.plot([1,1],[0,np.max(rho)],'k--')
plt.xlim([0,2])
plt.ylim([0,np.max(rho)])
plt.xlabel(r"Paramètre de relaxation $\omega$",fontsize=12)
plt.ylabel(r"Rayon spectral de la matrice d'itération $\rho(C)$",fontsize=12)
plt.title(r"Position ECEF du satellite GPS 1 = "+str(pos_sat1.astype(np.int64)),fontsize=12)
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap5_relaxation_optimale_s1_'+str(cas)+'.png')

#Méthode de la relaxation------------------------------------------------------

def relaxation(A,b,omega,x_0,n_max,e):
    
    #Récupérer les dimensions de la matrice A :
    m,n = np.shape(A)
    
    #Vérification des dimensions de A (nxn) et b (n) :
    if (m!=n)or(len(b)!=n):
        
        raise ValueError("Le système n'est pas de Cramer")
        
    #Récupérer les valeurs sur la diagonale de A :
    A_diag = np.diag(A)
    
    #Créer la matrice diagonale D, ayant les mêmes valeurs que la diagonale de A :
    D = np.diag(A_diag)
    
    #Créer la matrice A-D, ayant des zéros sur sa diagonale et les mêmes valeurs
    #que A partout ailleurs :
    AD = A-D
        
    #Vérifier si la matrice A est à diagonale strictement dominante :
    for i in range(n):
        
        if sum(abs(AD[:,i]))>=abs(A_diag[i]):
            
            print("Attention : la matrice A n'est pas à diagonale strictement dominante")
            break
        
    #Initialisation des variables :
    n = 0 #Nombre d'itérations
    x_n_old = np.copy(x_0) #Estimation de la solution à l'itération n-1
    x_n = np.copy(x_0) #Estimation de la solution à l'itération n
    for i in range(len(b)):
        x_n[i] = omega*(b[i]-np.dot(A[i,:i],x_n[:i])-np.dot(A[i,i+1:],x_n_old[i+1:]))/A[i,i] + (1-omega)*x_n_old[i]
    r_n = np.dot(A,x_n)-b #Résidu
            
    #Itérations de l'algorithme de la relaxation
	#tant qu'une des conditions d'arrêt n'est pas atteinte :
    while (n<n_max)and(np.linalg.norm(x_n-x_n_old,ord=2)>e)and(np.linalg.norm(r_n,ord=2)>e):
        
        #Mettre à jour l'estimation de la solution :
        x_n_old = np.copy(x_n) #Itération n
        for i in range(len(b)):
            x_n[i] = omega*(b[i]-np.dot(A[i,:i],x_n[:i])-np.dot(A[i,i+1:],x_n_old[i+1:]))/A[i,i] + (1-omega)*x_n_old[i] #Iteration n+1
        
        #Incrémenter le nombre d'itérations :
        n+=1
        
        #Renvoyer l'estimation de la solution du système et le résidu :
        r_n = np.dot(A,x_n)-b
                                                                    
    #Renvoyer l'estimation de la solution du système et le résidu :
    return x_n,r_n

#Application-------------------------------------------------------------------

omega = 1.25

x_0 = np.array([0,0,0],dtype=np.float64)

x_n,r_n = relaxation(A,b,omega,x_0,100,1e-3)