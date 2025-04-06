import matplotlib.pyplot as plt
import numpy as np

N0 = 8000 #Paramètre de forme (m3/mm)
R = 5 #Taux de pluie (mm/h)

def f(D):
    
    return N0*np.exp(-(4.1*R**-0.21)*D)*D**6

D = np.linspace(0,6,1000)


M = 5 #Nombre de sous-intervalles

#Methodes composite------------------------------------------------------------

def methode_composite(f,a,b,methode,M):
    
    #Découpage de l'intervalle [a,b] en M sous-intervalles avec un pas de 
    #(b-a)/M :
    x_i = [a+i*(b-a)/M for i in range(M+1)]
    
    #Initialisation de la somme des aires sous la courbe des différents 
    #sous-intervalles :
    somme = 0 
    
    #Boucle sur les sous-intervalles :
    for i in range(M):
              
        #Addition au compteur de l'aire sous la courbe pour ce sous-intervalle :
        somme += methode(f,x_i[i],x_i[i+1])
    
    #Renvoyer l'estimation de l'aire sous la courbe pour l'intervalle [a,b] :
    return somme

x = [1+i*(3-1)/M for i in range(M+1)]

#Rectangles composite----------------------------------------------------------

#A gauche :
    
def rectangles_gauche(f,a,b):
    
    return (b-a)*f(a)
    
plt.figure(1)
plt.plot(D,f(D),'r-')
for i in range(M):
    plt.fill([x[i],x[i],x[i+1],x[i+1]],[0,f(x[i]),f(x[i]),0],'pink')
    plt.plot([x[i],x[i],x[i+1],x[i+1]],[0,f(x[i]),f(x[i]),0],color='deeppink',linestyle='dotted')
for i in range(M):
    plt.scatter(x[i],f(x[i]),color='r',marker='o')
plt.plot([1,1],[0,1600],'k--')
plt.plot([3,3],[0,1600],'k--')
plt.xlim([0,6])
plt.ylim([0,1600])
plt.xlabel("D = taille des gouttes de pluie [mm]",fontsize=14)
plt.ylabel("f(D) [mm5/m3]",fontsize=14)
plt.title('Méthode des rectangles composite à gauche (M = '+str(M)+')')
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap4_rectangles_composite_1.png')

print('Rectangles composite à gauche = '+str(methode_composite(f,1,3,rectangles_gauche,M)))

#A droite :

def rectangles_droite(f,a,b):
    
    return (b-a)*f(b)
    
plt.figure(2)
plt.plot(D,f(D),'r-')
for i in range(M):
    plt.fill([x[i],x[i],x[i+1],x[i+1]],[0,f(x[i+1]),f(x[i+1]),0],'pink')
    plt.plot([x[i],x[i],x[i+1],x[i+1]],[0,f(x[i+1]),f(x[i+1]),0],color='deeppink',linestyle='dotted')
for i in range(M):
    plt.scatter(x[i+1],f(x[i+1]),color='r',marker='o')
plt.plot([1,1],[0,1600],'k--')
plt.plot([3,3],[0,1600],'k--')
plt.xlim([0,6])
plt.ylim([0,1600])
plt.xlabel("D = taille des gouttes de pluie [mm]",fontsize=14)
plt.ylabel("f(D) [mm5/m3]",fontsize=14)
plt.title('Méthode des rectangles composite à droite (M = '+str(M)+')')
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap4_rectangles_composite_2.png')

print('Rectangles composite à droite = '+str(methode_composite(f,1,3,rectangles_droite,M)))

#Point milieu :
    
def rectangles_milieu(f,a,b):
    
    return (b-a)*f((a+b)/2)

plt.figure(3)
plt.plot(D,f(D),'r-')
for i in range(M):
    plt.fill([x[i],x[i],x[i+1],x[i+1]],[0,f((x[i]+x[i+1])/2),f((x[i]+x[i+1])/2),0],'pink')
    plt.plot([x[i],x[i],x[i+1],x[i+1]],[0,f((x[i]+x[i+1])/2),f((x[i]+x[i+1])/2),0],color='deeppink',linestyle='dotted')
for i in range(M):
    plt.scatter((x[i]+x[i+1])/2,f((x[i]+x[i+1])/2),color='r',marker='o')
plt.plot([1,1],[0,1600],'k--')
plt.plot([3,3],[0,1600],'k--')
plt.xlim([0,6])
plt.ylim([0,1600])
plt.xlabel("D = taille des gouttes de pluie [mm]",fontsize=14)
plt.ylabel("f(D) [mm5/m3]",fontsize=14)
plt.title('Méthode des rectangles composite au point milieu (M = '+str(M)+')')
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap4_rectangles_composite_3.png')

print('Rectangles composite au point milieu = '+str(methode_composite(f,1,3,rectangles_milieu,M)))

#Trapezes composite------------------------------------------------------------

def trapezes(f,a,b):
    
    return (b-a)*(f(a)+f(b))/2

plt.figure(4)
plt.plot(D,f(D),'r-')

plt.plot([1,1],[0,1600],'k--')
plt.plot([3,3],[0,1600],'k--')
plt.xlim([0,6])
plt.ylim([0,1600])
plt.xlabel("D = taille des gouttes de pluie [mm]",fontsize=14)
plt.ylabel("f(D) [mm5/m3]",fontsize=14)
plt.title('Méthode des trapèzes composite (M = '+str(M)+')')
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap4_trapezes_composite.png')

print('Trapèzes composite = '+str(methode_composite(f,1,3,trapezes,M)))

#Simpson composite-------------------------------------------------------------

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

plt.plot([1,1],[0,1600],'k--')
plt.plot([3,3],[0,1600],'k--')
plt.xlim([0,6])
plt.ylim([0,1600])
plt.xlabel("D = taille des gouttes de pluie [mm]",fontsize=14)
plt.ylabel("f(D) [mm5/m3]",fontsize=14)
plt.title('Méthode de Simpson composite (M = '+str(M)+')')
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap4_simpson_composite.png')

print('Simpson composite = '+str(methode_composite(f,1,3,simpson,M)))