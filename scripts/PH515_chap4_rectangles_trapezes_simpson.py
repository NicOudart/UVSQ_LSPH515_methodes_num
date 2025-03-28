import matplotlib.pyplot as plt
import numpy as np
import scipy

N0 = 8000 #Paramètre de forme (m3/mm)
R = 5 #Taux de pluie (mm/h)

def f(D):
    
    return N0*np.exp(-(4.1*R**-0.21)*D)*D**6

D = np.linspace(0,6,1000)

#Rectangles--------------------------------------------------------------------

#A gauche :
plt.figure(1)
plt.plot(D,f(D),'r-')
plt.fill([1,1,3,3],[0,f(1),f(1),0],'pink')
plt.scatter(1,f(1),color='r',marker='o')
plt.xlim([0,6])
plt.ylim([0,1600])
plt.xlabel("D = taille des gouttes de pluie [mm]",fontsize=14)
plt.ylabel("f(D) [mm5/m3]",fontsize=14)
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap4_rectangles_1.png')

#A droite :
plt.figure(2)
plt.plot(D,f(D),'r-')
plt.fill([1,1,3,3],[0,f(3),f(3),0],'pink')
plt.scatter(3,f(3),color='r',marker='o')
plt.xlim([0,6])
plt.ylim([0,1600])
plt.xlabel("D = taille des gouttes de pluie [mm]",fontsize=14)
plt.ylabel("f(D) [mm5/m3]",fontsize=14)
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap4_rectangles_2.png')

#Point milieu :
plt.figure(3)
plt.plot(D,f(D),'r-')
plt.fill([1,1,3,3],[0,f(2),f(2),0],'pink')
plt.scatter(2,f(2),color='r',marker='o')
plt.xlim([0,6])
plt.ylim([0,1600])
plt.xlabel("D = taille des gouttes de pluie [mm]",fontsize=14)
plt.ylabel("f(D) [mm5/m3]",fontsize=14)
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap4_rectangles_3.png')

#Trapezes----------------------------------------------------------------------

plt.figure(4)
plt.plot(D,f(D),'r-')
plt.fill([1,1,3,3],[0,f(1),f(3),0],'pink')
plt.scatter([1,3],[f(1),f(3)],color='r',marker='o')
plt.xlim([0,6])
plt.ylim([0,1600])
plt.xlabel("D = taille des gouttes de pluie [mm]",fontsize=14)
plt.ylabel("f(D) [mm5/m3]",fontsize=14)
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap4_trapezes.png')

#Simpson-----------------------------------------------------------------------

plt.figure(5)
plt.plot(D,f(D),'r-')
plt.scatter([1,2,3],[f(1),f(2),f(3)],color='r',marker='o')
plt.xlim([0,6])
plt.ylim([0,1600])
plt.xlabel("D = taille des gouttes de pluie [mm]",fontsize=14)
plt.ylabel("f(D) [mm5/m3]",fontsize=14)
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap4_simpson.png')