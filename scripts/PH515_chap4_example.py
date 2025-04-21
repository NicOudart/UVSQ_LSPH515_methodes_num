import matplotlib.pyplot as plt
import numpy as np
import scipy

N0 = 8000 #Paramètre de forme (m3/mm)
R = 5 #Taux de pluie (mm/h)

def f(D):
    
    return N0*np.exp(-(4.1*R**-0.21)*D)*D**6

D = np.linspace(0,6,6000)

plt.figure(0)
plt.plot(D,f(D),'r-')
plt.fill([1]+[D[i] for i in range(1000,3000)]+[3],[0]+[f(D[i]) for i in range(1000,3000)]+[0],'pink')
plt.xlim([0,6])
plt.ylim([0,1600])
plt.xlabel("D = taille des gouttes de pluie [mm]",fontsize=14)
plt.ylabel("f(D) [mm5/m3]",fontsize=14)
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap4_example_0.png')

#Valeur théorique :
D = np.linspace(1,3,1000000)
Z=scipy.integrate.simpson(f(D),x=D)