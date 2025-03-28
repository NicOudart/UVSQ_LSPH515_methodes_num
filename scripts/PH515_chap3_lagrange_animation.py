import matplotlib.pyplot as plt
import numpy as np

l = 48.81094*np.pi/180 #Latitude de l'UFR des sciences
a = 23.438403*np.pi/180 #Latitude des tropiques

def f(d):
    
    d = d*2*np.pi/365.25 #Convertion du jour en angle
    
    return 48/(2*np.pi)*np.arccos(np.tan(l)*np.tan(np.arcsin(np.sin(a)*np.sin(d))))

x = np.array([30,60,90,120,150,180,240,270,300,330],dtype='float')

xp = np.linspace(0,365,366)
yp = np.zeros(366)
    
n = len(x)

for i in range(n):
    
    Li = np.ones(366)

    for j in range(n):
        
        if j!=i:
            Li = Li*(xp-x[j])/(x[i]-x[j])
            
    Li = Li*f(x[i])
    
    yp = yp + Li
            
    plt.figure(i)
    plt.scatter(x,f(x),color='k',marker='o')
    plt.plot(xp,Li,'g-')
    plt.scatter(np.delete(x,i),np.zeros(n-1),color='g',marker='o')
    plt.scatter(x[i],f(x[i]),color='g',marker='o')
    plt.xlim([0,365])
    plt.ylim([-2,17])
    plt.xlabel("x = jours depuis l'équinoxe de printemps",fontsize=14)
    plt.ylabel("f(x) = durée du jour en heures",fontsize=14)
    plt.title('f(x'+str(i)+')L'+str(i)+'(x)')
    plt.grid()
    plt.tight_layout()
    
    plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap3_lagrange_animation_'+str(i)+'.png')

plt.figure(n)
plt.scatter(x,f(x),color='k',marker='o')
plt.plot(xp,yp,'r-')
plt.xlim([0,365])
plt.ylim([-2,17])
plt.xlabel("x = jours depuis l'équinoxe de printemps",fontsize=14)
plt.ylabel("f(x) = durée du jour en heures",fontsize=14)
plt.title('p(x) = f(x0)L0(x) + f(x1)L1(x) + ... + f(x9)L9(x)')
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap3_lagrange_animation_'+str(n)+'.png')

plt.scatter(210,f(210),color='r',marker='o')
plt.plot([0,210],[f(210),f(210)],'r--')
plt.plot([210,210],[-2,f(210)],'r--')

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap3_lagrange_animation_'+str(n+1)+'.png')
