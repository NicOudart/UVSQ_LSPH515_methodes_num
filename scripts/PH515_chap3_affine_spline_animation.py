import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

l = 48.81094*np.pi/180 #Latitude de l'UFR des sciences
a = 23.438403*np.pi/180 #Latitude des tropiques

def f(d):
    
    d = d*2*np.pi/365.25 #Convertion du jour en angle
    
    return 48/(2*np.pi)*np.arccos(np.tan(l)*np.tan(np.arcsin(np.sin(a)*np.sin(d))))

x = np.array([30,60,90,120,150,180,240,270,300,330],dtype='float')
y = f(x)
    
plt.figure(1)
plt.scatter(x,y,color='k',marker='o')
plt.xlim([0,365])
plt.ylim([7,17])
plt.xlabel("x = jours depuis l'équinoxe de printemps",fontsize=14)
plt.ylabel("f(x) = durée du jour en heures",fontsize=14)
plt.title('p(x) obtenu par interpolation affine')
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap3_affine_0.png')

plt.plot(x,y,'r-')

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap3_affine_1.png')

xp = 210
yp = f(180)+((f(240)-f(180))/(240-180))*(210-180)

plt.scatter(xp,yp,color='r',marker='o')
plt.plot([0,xp],[yp,yp],'r--')
plt.plot([xp,xp],[7,yp],'r--')

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap3_affine_2.png')

cs = CubicSpline(x, y, bc_type='natural')

xp = np.linspace(0,365,366)

plt.figure(2)
plt.scatter(x,y,color='k',marker='o')
plt.xlim([0,365])
plt.ylim([7,17])
plt.xlabel("x = jours depuis l'équinoxe de printemps",fontsize=14)
plt.ylabel("f(x) = durée du jour en heures",fontsize=14)
plt.title('p(x) obtenu par interpolation spline cubique')
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap3_spline_0.png')

plt.plot(xp,cs(xp),'r-')

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap3_spline_1.png')

plt.scatter(210,cs(210),color='r',marker='o')
plt.plot([0,210],[cs(210),cs(210)],'r--')
plt.plot([210,210],[7,cs(210)],'r--')

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap3_spline_2.png')