import matplotlib.pyplot as plt
import numpy as np

l = 48.81094*np.pi/180 #Latitude de l'UFR des sciences
a = 23.438403*np.pi/180 #Latitude des tropiques

def f(d):
    
    d = d*2*np.pi/365.25 #Convertion du jour en angle
    
    return 48/(2*np.pi)*np.arccos(np.tan(l)*np.tan(np.arcsin(np.sin(a)*np.sin(d))))

x = np.array([30,60,90,120,150,180,240,270,300,330],dtype='float')
y = f(x)

xp = np.linspace(0,365,366)
yp = np.zeros(366)
    
n = len(x)
    
c = np.zeros(n)

for i in range(n):
    
    c[i] = y[i]

for i in range(1,n):
    
    for k in range(n-1,i-1,-1):
        
        c[k] = (c[k]-c[k-1])/(x[k]-x[k-i])
        
yp = np.ones(366)*c[0]

titre = 'p(x) = c0'
        
plt.figure(0,figsize=(10, 5))
plt.scatter(x,f(x),color='k',marker='o')
plt.plot(xp,yp,'g-')
plt.scatter(x[0],f(x[0]),color='g',marker='o')
plt.xlim([0,365])
plt.ylim([-2,17])
plt.xlabel("x = jours depuis l'équinoxe d'automne",fontsize=14)
plt.ylabel("f(x) = durée du jour en heures",fontsize=14)
plt.title(titre)
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap3_Newton_animation_0.png')  

p = c[0]     

for i in range(1,n):
    
    titre += ' + c'+str(i)+'v'+str(i)+'(x)'
    
    v = c[i]
    
    for j in range(i):
        
        v = v*(xp-x[j])
        
    p += v
    
    plt.figure(i,figsize=(10, 5))
    plt.scatter(x,f(x),color='k',marker='o')
    plt.plot(xp,p,'g-')
    plt.scatter(x[:i+1],f(x[:i+1]),color='g',marker='o')
    plt.xlim([0,365])
    plt.ylim([-2,17])
    plt.xlabel("x = jours depuis l'équinoxe d'automne",fontsize=14)
    plt.ylabel("f(x) = durée du jour en heures",fontsize=14)
    plt.title(titre)
    plt.grid()
    plt.tight_layout()
    
    plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap3_Newton_animation_'+str(i)+'.png')
    
plt.figure(n,figsize=(10, 5))
plt.scatter(x,f(x),color='k',marker='o')
plt.plot(xp,p,'r-')
plt.xlim([0,365])
plt.ylim([-2,17])
plt.xlabel("x = jours depuis l'équinoxe d'automne",fontsize=14)
plt.ylabel("f(x) = durée du jour en heures",fontsize=14)
plt.title(titre)
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap3_Newton_animation_'+str(n)+'.png')

plt.scatter(210,f(210),color='r',marker='o')
plt.plot([0,210],[f(210),f(210)],'r--')
plt.plot([210,210],[-2,f(210)],'r--')

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap3_Newton_animation_'+str(n+1)+'.png')
