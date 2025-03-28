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
plt.plot([1,1],[0,1600],'k--')
plt.plot([3,3],[0,1600],'k--')
plt.scatter(1,f(1),color='r',marker='o')
plt.annotate('a',(1,1500),textcoords="offset points",xytext=(-10,0),ha='center',fontsize=14,color='k')
plt.annotate('b',(3,1500),textcoords="offset points",xytext=(10,0),ha='center',fontsize=14,color='k')
plt.annotate('f(x_0)',(1,f(1)),textcoords="offset points",xytext=(-30,10),ha='center',fontsize=14,color='r')
plt.xlim([0,6])
plt.ylim([0,1600])
plt.xlabel("D = taille des gouttes de pluie [mm]",fontsize=14)
plt.ylabel("f(D) [mm5/m3]",fontsize=14)
plt.title('Méthode des rectangles à gauche')
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap4_rectangles_1.png')

#A droite :
plt.figure(2)
plt.plot(D,f(D),'r-')
plt.fill([1,1,3,3],[0,f(3),f(3),0],'pink')
plt.plot([1,1],[0,1600],'k--')
plt.plot([3,3],[0,1600],'k--')
plt.scatter(3,f(3),color='r',marker='o')
plt.annotate('a',(1,1500),textcoords="offset points",xytext=(-10,0),ha='center',fontsize=14,color='k')
plt.annotate('b',(3,1500),textcoords="offset points",xytext=(10,0),ha='center',fontsize=14,color='k')
plt.annotate('f(x_0)',(3,f(3)),textcoords="offset points",xytext=(30,10),ha='center',fontsize=14,color='r')
plt.xlim([0,6])
plt.ylim([0,1600])
plt.xlabel("D = taille des gouttes de pluie [mm]",fontsize=14)
plt.ylabel("f(D) [mm5/m3]",fontsize=14)
plt.title('Méthode des rectangles à droite')
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap4_rectangles_2.png')

#Point milieu :
plt.figure(3)
plt.plot(D,f(D),'r-')
plt.fill([1,1,3,3],[0,f(2),f(2),0],'pink')
plt.plot([1,1],[0,1600],'k--')
plt.plot([3,3],[0,1600],'k--')
plt.scatter(2,f(2),color='r',marker='o')
plt.annotate('a',(1,1500),textcoords="offset points",xytext=(-10,0),ha='center',fontsize=14,color='k')
plt.annotate('b',(3,1500),textcoords="offset points",xytext=(10,0),ha='center',fontsize=14,color='k')
plt.annotate('f(x_0)',(2,f(2)),textcoords="offset points",xytext=(5,-30),ha='center',fontsize=14,color='r')
plt.xlim([0,6])
plt.ylim([0,1600])
plt.xlabel("D = taille des gouttes de pluie [mm]",fontsize=14)
plt.ylabel("f(D) [mm5/m3]",fontsize=14)
plt.title('Méthode des rectangles au point milieu')
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap4_rectangles_3.png')

#Trapezes----------------------------------------------------------------------

plt.figure(4)
plt.plot(D,f(D),'r-')
plt.fill([1,1,3,3],[0,f(1),f(3),0],'pink')
plt.plot([1,1],[0,1600],'k--')
plt.plot([3,3],[0,1600],'k--')
plt.scatter([1,3],[f(1),f(3)],color='r',marker='o')
plt.annotate('a',(1,1500),textcoords="offset points",xytext=(-10,0),ha='center',fontsize=14,color='k')
plt.annotate('b',(3,1500),textcoords="offset points",xytext=(10,0),ha='center',fontsize=14,color='k')
plt.annotate('f(a)',(1,f(1)),textcoords="offset points",xytext=(-20,10),ha='center',fontsize=14,color='r')
plt.annotate('f(b)',(3,f(3)),textcoords="offset points",xytext=(20,10),ha='center',fontsize=14,color='r')
plt.xlim([0,6])
plt.ylim([0,1600])
plt.xlabel("D = taille des gouttes de pluie [mm]",fontsize=14)
plt.ylabel("f(D) [mm5/m3]",fontsize=14)
plt.title('Méthode des trapèzes')
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap4_trapezes.png')

#Simpson-----------------------------------------------------------------------

plt.figure(5)
plt.plot(D,f(D),'r-')
plt.plot([1,1],[0,1600],'k--')
plt.plot([3,3],[0,1600],'k--')
plt.scatter([1,2,3],[f(1),f(2),f(3)],color='r',marker='o')
plt.annotate('a',(1,1500),textcoords="offset points",xytext=(-10,0),ha='center',fontsize=14,color='k')
plt.annotate('b',(3,1500),textcoords="offset points",xytext=(10,0),ha='center',fontsize=14,color='k')
plt.annotate('f(a)',(1,f(1)),textcoords="offset points",xytext=(-20,10),ha='center',fontsize=14,color='r')
plt.annotate('f(b)',(3,f(3)),textcoords="offset points",xytext=(20,10),ha='center',fontsize=14,color='r')
plt.xlim([0,6])
plt.ylim([0,1600])
plt.xlabel("D = taille des gouttes de pluie [mm]",fontsize=14)
plt.ylabel("f(D) [mm5/m3]",fontsize=14)
plt.title('Méthode de Simpson')
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap4_simpson.png')