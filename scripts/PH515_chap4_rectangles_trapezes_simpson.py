import matplotlib.pyplot as plt
import numpy as np

N0 = 8000 #Paramètre de forme (m3/mm)
R = 5 #Taux de pluie (mm/h)

def f(D):
    
    return N0*np.exp(-(4.1*R**-0.21)*D)*D**6

D = np.linspace(0,6,1000)

#Rectangles--------------------------------------------------------------------

#A gauche :
    
def rectangles_gauche(f,a,b):
    
    return (b-a)*f(a)
    
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

print('Rectangles à gauche = '+str(rectangles_gauche(f,1,3)))

#A droite :

def rectangles_droite(f,a,b):
    
    return (b-a)*f(b)
    
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

print('Rectangles à droite = '+str(rectangles_droite(f,1,3)))

#Point milieu :
    
def rectangles_milieu(f,a,b):
    
    return (b-a)*f((a+b)/2)

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

print('Rectangles au point milieu = '+str(rectangles_milieu(f,1,3)))

#Trapezes----------------------------------------------------------------------

def trapezes(f,a,b):
    
    return (b-a)*(f(a)+f(b))/2

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

print('Trapèzes = '+str(trapezes(f,1,3)))

#Simpson-----------------------------------------------------------------------

def simpson(f,a,b):
    
    return (b-a)*(f(a)+4*f((a+b)/2)+f(b))/6

#Interpolation de Lagrange pour déterminer la parabole :
def lagrange(x,y,xp):
    n = len(x)
    yp = 0
    for i in range(n):
        Li = 1
        for j in range(n):
            if j!=i:
                Li = Li*(xp-x[j])/(x[i]-x[j])
        yp = yp + y[i]*Li
    return yp

D_fill = np.linspace(1,3,3000)

plt.figure(5)
plt.plot(D,f(D),'r-')
plt.fill([1]+[d for d in D_fill]+[3],[0]+[lagrange([1,2,3],[f(1),f(2),f(3)],d) for d in D_fill]+[0],'pink')
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

print('Simpson = '+str(simpson(f,1,3)))