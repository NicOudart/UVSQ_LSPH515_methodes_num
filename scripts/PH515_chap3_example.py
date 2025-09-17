import matplotlib.pyplot as plt
import numpy as np

l = 48.81094*np.pi/180 #Latitude de l'UFR des sciences
a = 23.438403*np.pi/180 #Latitude des tropiques

def f(d):
    
    return 48/(2*np.pi)*np.arccos(np.tan(l)*np.tan(np.arcsin(np.sin(a)*np.sin(d*2*np.pi/365.25))))

days = np.array([30,60,90,120,150,180,240,270,300,330],dtype='float')

plt.figure(0)
plt.scatter(days,f(days),color='k',marker='o')
plt.xlim([0,365])
plt.ylim([7,17])
plt.xlabel("x = jours depuis l'équinoxe d'automne",fontsize=14)
plt.ylabel("f(x) = durée du jour en heures",fontsize=14)
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap3_example_animation_0.png')

plt.scatter(210,f(210),color='r',marker='o')

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap3_example_animation_1.png')