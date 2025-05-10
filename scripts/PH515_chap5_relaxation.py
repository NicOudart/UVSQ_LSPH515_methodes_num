import numpy as np
import matplotlib.pyplot as plt

xs_1 = 22250

pos_rec = np.array([4205,158,4777],dtype=np.float64) #Coordonnées du récepteur GPS

pos_sat1 = np.array([xs_1,13000,18000],dtype=np.float64) #Coordonnées du satellite GPS 1
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
E = np.tril(A)-D
F = np.triu(A)-D

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
plt.xlabel(r"Paramètre de relaxation $\omega$",fontsize=13)
plt.ylabel(r"Rayon spectral de la matrice d'itération $\rho(C)$",fontsize=12)
plt.title(r"Recherche du paramètre de relaxation optimal ($x_{s1}$ = "+str(xs_1)+")",fontsize=13)
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap5_relaxation_optimale_xs1_'+str(xs_1)+'.png')
