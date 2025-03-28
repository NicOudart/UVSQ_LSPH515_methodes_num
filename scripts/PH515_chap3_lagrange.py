import numpy as np

l = 48.81094*np.pi/180 #Latitude de l'UFR des sciences
a = 23.438403*np.pi/180 #Latitude des tropiques

def f(d):
    
    d = d*2*np.pi/365.25 #Convertion du jour en angle
    
    return 48/(2*np.pi)*np.arccos(np.tan(l)*np.tan(np.arcsin(np.sin(a)*np.sin(d))))

x = np.array([30,60,90,120,150,180,240,270,300,330],dtype='float')

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

yp = lagrange(x,f(x),210)
    
