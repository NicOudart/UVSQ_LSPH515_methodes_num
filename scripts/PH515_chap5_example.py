import numpy as np
import scipy

pos_rec = np.array([4205,158,4777]) #Coordonnées du récepteur GPS

pos_sat1 = np.array([14000,3000,25000]) #Coordonnées du satellite GPS 1
pos_sat2 = np.array([24000,-3000,15000]) #Coordonnées du satellite GPS 2
pos_sat3 = np.array([18000,-8000,21000]) #Coordonnées du satellite GPS 3
pos_sat4 = np.array([22000,14000,12000]) #Coordonnées du satellite GPS 4

#Génération des temps de retard------------------------------------------------

t1 = ((sum((pos_rec-pos_sat1)**2))**0.5)/3e5
t2 = ((sum((pos_rec-pos_sat2)**2))**0.5)/3e5
t3 = ((sum((pos_rec-pos_sat3)**2))**0.5)/3e5
t4 = ((sum((pos_rec-pos_sat4)**2))**0.5)/3e5

#Génération des matrices A et b------------------------------------------------

A = np.vstack((pos_sat2-pos_sat1,pos_sat3-pos_sat1,pos_sat4-pos_sat1))

b_row1 = 0.5*(sum(pos_sat2**2)-sum(pos_sat1**2)-(3e5*t2)**2+(3e5*t1)**2)
b_row2 = 0.5*(sum(pos_sat3**2)-sum(pos_sat1**2)-(3e5*t3)**2+(3e5*t1)**2)
b_row3 = 0.5*(sum(pos_sat4**2)-sum(pos_sat1**2)-(3e5*t4)**2+(3e5*t1)**2)

b = np.array([b_row1,b_row2,b_row3])

#Solution théorique :
pos_mes = scipy.linalg.solve(A,b)