import matplotlib.pyplot as plt
import numpy as np

#Exemple de déplacement du ballon----------------------------------------------

#Fonction du déplacement en fonction de t :
def p(t):
    
    return 100/(1+np.exp(-0.1*t+5))

#Fonction de la vitesse théorique en fonction de t :
def w(t):
    
    p_t = p(t)
    
    return 0.1*p_t*(100-p_t)/100

#Methode de la différence décentrée à droite-----------------------------------

#Définition de la méthode sous la forme d'une fonction :
def difference_decentree_droite(f,h):

	#On s'assure que f contient des réels encodés sur 64 bits :
	f.astype(dtype=np.float64)

	#On récupère le nombre d'éléments N de f :
	N = len(f)
	
	#On initialise un vecteur df de taille N-1 ne contenant que des zéros,
	#qui contiendra les valeurs de la dérivée de f :
	df = np.zeros(N-1,dtype=np.float64)
	
	#On fait une boucle sur les différentes valeurs de dérivée à calculer :
	for i in range(N-1):
	
		#On approxime la dérivée en un point par la formule de la différence décentrée à droite :
		df[i] = (f[i+1]-f[i])/h
	
	#Renvoyer le vecteur des valeurs de la dérivée de f :
	return df

#Application de la méthode à l'exemple-----------------------------------------

#Evaluations dans le cas "continu" pour référence :
t_c = np.linspace(0,100,1000) #Axe du temps
p_c = p(t_c) #Déplacement théorique du ballon
w_c = w(t_c) #Vitesse théorique du vent

#Evaluations dans le cas discret :
h = 10 #Pas de discrétisation
t_i = np.arange(0,101,h) #Instants pour lesquels le déplacement du ballon est mesuré
p_i = p(t_i) #Déplacement mesuré du ballon
w_i = difference_decentree_droite(p_i,h) #Estimation de la vitesse du vent

#Affichage graphique-----------------------------------------------------------

#Problème initial :
    
plt.figure(0,figsize=(10, 5))

plt.subplot(1, 2, 1)
plt.plot(t_c,p_c,'k--')
plt.scatter(t_i,p_i,color='k',marker='o')
plt.xlim([0,100])
plt.ylim([0,110])
plt.xlabel("t [s]",fontsize=14)
plt.ylabel("p(t) [m]",fontsize=14)
plt.title('Déplacement du ballon')
plt.grid()
plt.tight_layout()

plt.subplot(1, 2, 2)
plt.plot(t_c,w_c,'k--')
plt.xlim([0,100])
plt.ylim([0,2.75])
plt.xlabel("t [s]",fontsize=14)
plt.ylabel("W(t) [m/s]",fontsize=14)
plt.title('Vitesse du vent')
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap1_exemple_0.png')

#Application de la méthode :

N = len(t_i)

for i in range(N-1):

    plt.figure(i+1,figsize=(10, 5))

    plt.subplot(1, 2, 1)
    plt.plot(t_c,p_c,'k--')
    plt.scatter(t_i,p_i,color='k',marker='o')
    plt.plot(t_c,w_i[i]*(t_c-t_i[i])+p_i[i],'r--')
    plt.scatter(t_i[i:i+2],p_i[i:i+2],color='r',marker='o')
    plt.xlim([0,100])
    plt.ylim([0,110])
    plt.xlabel("t [s]",fontsize=14)
    plt.ylabel("p(t) [m]",fontsize=14)
    plt.title('Déplacement du ballon - p(t_'+str(i)+')')
    plt.grid()
    plt.tight_layout()

    plt.subplot(1, 2, 2)
    plt.plot(t_c,w_c,'k--')
    plt.scatter(t_i[:i+1],w_i[:i+1],color='r',marker='o')
    plt.xlim([0,100])
    plt.ylim([0,2.75])
    plt.xlabel("t [s]",fontsize=14)
    plt.ylabel("W(t) [m/s]",fontsize=14)
    plt.title('Vitesse du vent - W(t_'+str(i)+')')
    plt.grid()
    plt.tight_layout()
    
    plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap1_exemple_'+str(i+1)+'.png')

#Solution du problème et erreurs :

plt.figure(N,figsize=(10, 5))

plt.subplot(1, 2, 1)
plt.plot(t_c,p_c,'k--')
plt.scatter(t_i,p_i,color='k',marker='o')
plt.xlim([0,100])
plt.ylim([0,110])
plt.xlabel("t [s]",fontsize=14)
plt.ylabel("p(t) [m]",fontsize=14)
plt.title('Déplacement du ballon')
plt.grid()
plt.tight_layout()

plt.subplot(1, 2, 2)
plt.plot(t_c,w_c,'k--')
plt.scatter(t_i[:-1],w_i,color='r',marker='o')
for i in range(N-1):
    plt.plot([t_i[i],t_i[i]],[w_i[i],w(t_i[i])],'r--')
plt.xlim([0,100])
plt.ylim([0,2.75])
plt.xlabel("t [s]",fontsize=14)
plt.ylabel("W(t) [m/s]",fontsize=14)
plt.title('Vitesse du vent')
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap1_exemple_'+str(N)+'.png')