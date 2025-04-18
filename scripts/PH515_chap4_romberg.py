import matplotlib.pyplot as plt
import numpy as np

N0 = 8000 #Paramètre de forme (m3/mm)
R = 5 #Taux de pluie (mm/h)

def f(D):
    
    return N0*np.exp(-(4.1*R**-0.21)*D)*D**6

D = np.linspace(0,6,1000)

#Methode des trapèzes composite------------------------------------------------

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

#Accélération de Romberg:------------------------------------------------------

def romberg_trapezes(f,a,b,epsilon):
    
    #Initialiser le nombre de sous-intervalles d'intégration M :
    M = 1
    
    #Initialiser les aires estimées pour M et 2M sous-intervalles d'intégration :
    aire_M = trapezes_composite(f,a,b,M)
    aire_2M = trapezes_composite(f,a,b,2*M)
    
    #Boucler tant que la condition sur l'erreur d'intégration n'est pas atteinte :
    while abs(aire_M-aire_2M)>(3*epsilon):
        
        #Multiplier le nombre de sous-intervalles d'intégration par 2 :
        M *= 2
        
        #Estimer l'aire pour M et 2M sous-intervalles d'intégration :
        aire_M = trapezes_composite(f,a,b,M)
        aire_2M = trapezes_composite(f,a,b,2*M)
        
    return aire_2M,2*M

#Application au problème-------------------------------------------------------

epsilon = 1e-3 #Précision choisie

aire_M,M = romberg_trapezes(f,1,3,epsilon)

print('Intégrale estimée = '+str(aire_M)+', Nombre de sous-intervalles = '+str(M))

#Affichage---------------------------------------------------------------------

aire_ref = trapezes_composite(f,1,3,32768)

for idx_fig in range(6):
    
    nb_intervalles = 2**idx_fig
    
    x = [1+i*(3-1)/nb_intervalles for i in range(nb_intervalles+1)]
    
    aire_nb_intervalles = trapezes_composite(f,1,3,nb_intervalles)
    
    plt.figure(idx_fig)
    plt.plot(D,f(D),'r-')
    for i in range(nb_intervalles):
        plt.fill([x[i],x[i],x[i+1],x[i+1]],[0,f(x[i]),f(x[i+1]),0],'pink')
        plt.plot([x[i],x[i],x[i+1],x[i+1]],[0,f(x[i]),f(x[i+1]),0],color='deeppink',linestyle='dotted')
    for i in range(nb_intervalles):
        plt.scatter(x[i],f(x[i]),color='r',marker='o')
        plt.scatter(x[i+1],f(x[i+1]),color='r',marker='o')
    plt.plot([1,1],[0,1600],'k--')
    plt.plot([3,3],[0,1600],'k--')
    plt.xlim([0,6])
    plt.ylim([0,1600])
    plt.xlabel("D = taille des gouttes de pluie [mm]",fontsize=14)
    plt.ylabel("f(D) [mm5/m3]",fontsize=14)
    plt.title('Méthode des trapèzes composite (M = '+str(nb_intervalles)+')')
    plt.grid()
    plt.text(3.5,1300,"$\epsilon$ = "+str(round(abs(aire_nb_intervalles-aire_ref),2))+" mm6/m3",fontsize=12,backgroundcolor='pink')
    plt.tight_layout()
    
    plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap4_romberg_trapezes_'+str(idx_fig)+'.png')