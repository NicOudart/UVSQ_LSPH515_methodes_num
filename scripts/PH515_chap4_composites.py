import matplotlib.pyplot as plt
import numpy as np

N0 = 8000 #Paramètre de forme (m3/mm)
R = 5 #Taux de pluie (mm/h)

def f(D):
    
    return N0*np.exp(-(4.1*R**-0.21)*D)*D**6

D = np.linspace(0,6,1000)

M = 5 #Nombre de sous-intervalles

#Methodes simples--------------------------------------------------------------

def rectangles_gauche(f,a,b):
    
    return (b-a)*f(a)

def rectangles_droite(f,a,b):
    
    return (b-a)*f(b)

def rectangles_milieu(f,a,b):
    
    return (b-a)*f((a+b)/2)

def trapezes(f,a,b):
    
    return (b-a)*(f(a)+f(b))/2

def simpson(f,a,b):
    
    return (b-a)*(f(a)+4*f((a+b)/2)+f(b))/6

#Methodes composite------------------------------------------------------------

def methode_composite(f,a,b,methode,M):
    
    #Découpage de l'intervalle [a,b] en M sous-intervalles avec un pas de 
    #(b-a)/M :
    x_i = [a+i*(b-a)/M for i in range(M+1)]
    
    #Initialisation de la somme des aires sous la courbe des différents 
    #sous-intervalles :
    aire = 0 
    
    #Boucle sur les sous-intervalles :
    for i in range(M):
              
        #Addition au compteur de l'aire sous la courbe pour ce sous-intervalle :
        aire += methode(f,x_i[i],x_i[i+1])
    
    #Renvoyer l'estimation de l'aire sous la courbe pour l'intervalle [a,b] :
    return aire

def rectangles_gauche_composite(f,a,b,M):
    
    #Déterminer la largeur h de chaque sous-intervalle :
    h = (b-a)/M
    
    #Découpage de l'intervalle [a,b] en M sous-intervalles avec un pas de h :
    x_i = [a+i*h for i in range(M+1)]
    
    #Initialiser la somme des évaluations de f pour chaque sous-intervalle:
    somme = 0 
    
    #Boucle sur les sous-intervalles :
    for i in range(M):
        
        #Sommer la valeur de f "à gauche" du sous-intervalle :
        somme += f(x_i[i])
        
    #Calcul de la somme des aires des sous-intervalles :
    aire = somme*h
    
    return aire

def rectangles_droite_composite(f,a,b,M):
    
    #Déterminer la largeur h de chaque sous-intervalle :
    h = (b-a)/M
    
    #Découpage de l'intervalle [a,b] en M sous-intervalles avec un pas de h :
    x_i = [a+i*h for i in range(M+1)]
    
    #Initialiser la somme des évaluations de f pour chaque sous-intervalle:
    somme = 0 
    
    #Boucle sur les sous-intervalles :
    for i in range(M):
        
        #Sommer la valeur de f "à droite" du sous-intervalle :
        somme += f(x_i[i+1])
        
    #Calcul de la somme des aires des sous-intervalles :
    aire = somme*h
    
    return aire

def rectangles_milieu_composite(f,a,b,M):
    
    #Déterminer la largeur h de chaque sous-intervalle :
    h = (b-a)/M
    
    #Découpage de l'intervalle [a,b] en M sous-intervalles avec un pas de h :
    x_i = [a+i*h for i in range(M+1)]
    
    #Initialiser la somme des évaluations de f pour chaque sous-intervalle:
    somme = 0 
    
    #Boucle sur les sous-intervalles :
    for i in range(M):
        
        #Sommer la valeur de f "au milieu" du sous-intervalle :
        somme += f((x_i[i]+x_i[i+1])/2)
        
    #Calcul de la somme des aires des sous-intervalles :
    aire = somme*h
    
    return aire

def trapezes_composite(f,a,b,M):
    
    #Déterminer la largeur h de chaque sous-intervalle :
    h = (b-a)/M
    
    #Découpage de l'intervalle [a,b] en M sous-intervalles avec un pas de h :
    x_i = [a+i*h for i in range(M+1)]
    
    #Initialiser la somme des évaluations de f pour chaque sous-intervalle avec
    #f(a) et f(b):
    somme = (f(x_i[0])+f(x_i[-1]))/2
    
    #Boucle sur les sous-intervalles :
    for i in range(1,M):
        
        #Sommer la valeur de f à gauche du sous-intervalle :
        somme += f(x_i[i])
        
    #Calcul de la somme des aires des sous-intervalles :
    aire = somme*h
    
    return aire

def simpson_composite(f,a,b,M):
    
    #Déterminer la largeur h de chaque sous-intervalle :
    h = (b-a)/(2*M)
    
    #Découpage de l'intervalle [a,b] en M sous-intervalles avec un pas de h :
    x_i = [a+i*h for i in range(2*M+1)]
    
    #Initialiser la somme des évaluations de f pour les éléments pairs et 
    #impairs de l'intervalle :
    somme_pair = 0
    somme_impair = 0
    
    #1ère boucle sur les éléments pairs des sous-intervalles :
    for i in range(1,M):
        
        #Sommer la valeur de f du sous-intervalle :
        somme_pair += f(x_i[2*i])
    
    #2nde boucle sur les éléments impairs des sous-intervalles :
    for i in range(M):
        somme_impair += f(x_i[2*i+1])
       
    #Additionner les valeurs de f :
    somme = f(x_i[0])+f(x_i[-1])+2*somme_pair+4*somme_impair
        
    #Calcul de la somme des aires des sous-intervalles :
    aire = somme*h/3
    
    return aire

#Vecteur de valeurs à évaluer pour l'affichage---------------------------------

x = [1+i*(3-1)/M for i in range(M+1)]

#Rectangles composite----------------------------------------------------------

#A gauche :
    
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

plt.figure(4)
plt.plot(D,f(D),'r-')
for i in range(M):
    plt.fill([x[i],x[i],x[i+1],x[i+1]],[0,f(x[i]),f(x[i+1]),0],'pink')
    plt.plot([x[i],x[i],x[i+1],x[i+1]],[0,f(x[i]),f(x[i+1]),0],color='deeppink',linestyle='dotted')
for i in range(M):
    plt.scatter(x[i],f(x[i]),color='r',marker='o')
    plt.scatter(x[i+1],f(x[i+1]),color='r',marker='o')
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
for i in range(M):
    D_fill = np.linspace(x[i],x[i+1],100)
    plt.fill([x[i]]+[d for d in D_fill]+[x[i+1]],[0]+[lagrange([x[i],(x[i]+x[i+1])/2,x[i+1]],[f(x[i]),f((x[i]+x[i+1])/2),f(x[i+1])],d) for d in D_fill]+[0],'pink')
    plt.plot([x[i]]+[d for d in D_fill]+[x[i+1]],[0]+[lagrange([x[i],(x[i]+x[i+1])/2,x[i+1]],[f(x[i]),f((x[i]+x[i+1])/2),f(x[i+1])],d) for d in D_fill]+[0],color='deeppink',linestyle='dotted')
for i in range(M):
    plt.scatter(x[i],f(x[i]),color='r',marker='o')
    plt.scatter((x[i]+x[i+1])/2,f((x[i]+x[i+1])/2),color='r',marker='o')
    plt.scatter(x[i+1],f(x[i+1]),color='r',marker='o')
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