import numpy as np

l = 48.81094*np.pi/180 #Latitude de l'UFR des sciences
a = 23.438403*np.pi/180 #Latitude des tropiques

def f(d):
    
    d = d*2*np.pi/365.25 #Convertion du jour en angle
    
    return 48/(2*np.pi)*np.arccos(np.tan(l)*np.tan(np.arcsin(np.sin(a)*np.sin(d))))

x = np.array([30,60,90,120,150,180,240,270,300,330],dtype='float')

def newton(x,y,xp):
    
    #Récupération du nombre de points connus :
    n = len(x)
    
    #Initialisation d'un vecteur nul qui contiendra les coefficients du 
    #polynôme de Newton :
    c = np.zeros(n)
    
    #Boucle sur les ordonnées connues pour déterminer la 1ère colonne du 
    #tableau des différences divisées (ordre 0) :
    for i in range(n):
        
        c[i] = y[i]
    
    #Boucle pour calculer les colonnes du tableau des différences divisées 
    #(de gauche à droite):
    for i in range(1,n):
        
        #Boucle pour calculer les lignes du tableau des différences divisées 
        #(du bas jusqu'à la diagonale de chaque colonnes):
        for k in range(n-1,i-1,-1):
            
            c[k] = (c[k]-c[k-1])/(x[k]-x[k-i])
            #(Avec cette formule récursive sur les éléments du vecteur des 
            #coefficients, on ne gardera en mémoire que les éléments de la 
            #diagonale du tableau des différences divisées).
          
    
    #Calcul de l'ordonnée du point interpolé avec l'algorithme de Horner :
    yp = c[n-1]
    
    for i in range(n-2,-1,-1):
        
        yp = c[i] + (xp-x[i])*yp
        
    #Renvoyer l'ordonnée du point interpolé :    
    return yp

yp = newton(x,f(x),210)
    
