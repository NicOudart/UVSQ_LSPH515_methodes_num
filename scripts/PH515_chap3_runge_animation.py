import matplotlib.pyplot as plt
import numpy as np

nb_interp = 50

def f(x):
    
    return 1/(1+32*(x**2))

def lagrange(x,y,xp):
    
    #Récupération du nombre de points connus :
    n = len(x)
    
    #Initialisation de l'ordonnée du point interpolé yp :
    yp = 0
    
    #1ère boucle sur les points connus (x[i],y[i]) :
    for i in range(n):
        
        #Initialisation du polynôme de la base de Lagrange associé au i-ème 
        #point connu Li :
        Li = 1

        #2ème boucle sur les abscisses connues x[j]:
        for j in range(n):
            
            #Dans le cas où les 2 points connus x[i] et x[j] sont distincts,
            #multiplication de Li par un coefficient obtenu avec xp, x[i] et 
            #x[j], de telle façon que Li = 1 si xp = x[i] et Li = 0 si xp = x[j] 
            #(voir formule du cours) :
            if j!=i:
                Li = Li*(xp-x[j])/(x[i]-x[j])
           
        #Addition à yp de la valeur du polynôme y[i]xLi, qui est égal à y[i] en 
        #x[i] et nul pour tout les x[j] (voir formule du cours):
        yp = yp + y[i]*Li
        
    #Renvoyer l'ordonnée du point interpolé :    
    return yp

x = np.linspace(-1,1,nb_interp)
y = f(x)

xp = np.linspace(-1,1,1000)
yp = lagrange(x,y,xp)

plt.figure()
plt.scatter(x,f(x),color='k',marker='o')
plt.plot(xp,yp,'r-')
plt.xlim([-1.05,1.05])
plt.ylim([-1.05,1.05])
plt.xlabel("x",fontsize=14)
plt.ylabel("f(x)",fontsize=14)
plt.title("p(x) avec "+str(nb_interp)+" points d'interpolation")
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap3_runge_animation_'+str(nb_interp)+'.png')

plt.figure()
plt.plot(xp,f(xp),'k-')
plt.xlim([-1.05,1.05])
plt.ylim([-1.05,1.05])
plt.xlabel("x",fontsize=14)
plt.ylabel("f(x)",fontsize=14)
plt.title("f(x) à interpoler")
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap3_runge_animation_0.png')

