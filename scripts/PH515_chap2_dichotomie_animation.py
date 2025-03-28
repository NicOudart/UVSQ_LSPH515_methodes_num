import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return (x**2)-2

x = np.linspace(-0.5,2.5,1000)

a_n = 1
b_n = 2 
x_n = (a_n+b_n)/2

#Step 0:----------------------------------

plt.figure(0)
plt.plot(x,f(x),'r-')
plt.plot([-0.5,2**0.5],[0,0],'k--')
plt.plot([2**0.5,2**0.5],[-3,0],'k--')
plt.scatter(2**0.5,0,color='k',marker='o')
plt.xlim([-0.5,2.5])
plt.ylim([-3,3])
plt.xlabel('x',fontsize=14)
plt.ylabel('f(x)',fontsize=14)
plt.grid()
plt.tight_layout()

plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap2_animation_0.png')

#Loop:----------------------------------

for idx in range(1,6):

    plt.figure(idx)
    plt.plot(x,f(x),'r-')
    plt.scatter(2**0.5,0,color='r',marker='o')
    plt.scatter(x_n,f(x_n),color='g',marker='o')
    plt.axvspan(-0.5,a_n,facecolor='k',alpha=0.2)
    plt.axvspan(b_n,2.5,facecolor='k',alpha=0.2)
    plt.scatter([a_n,b_n],[f(a_n),f(b_n)],color='k',marker='o')
    plt.annotate('x_'+str(idx-1),(x_n,f(x_n)),textcoords="offset points",xytext=(-5,10),ha='center',fontsize=14,color='g')
    plt.annotate('a_'+str(idx-1),(a_n,f(a_n)),textcoords="offset points",xytext=(-25,0),ha='center',fontsize=14)
    plt.annotate('b_'+str(idx-1),(b_n,f(b_n)),textcoords="offset points",xytext=(25,0),ha='center',fontsize=14)
    plt.xlim([-0.5,2.5])
    plt.ylim([-3,3])
    plt.xlabel('x',fontsize=14)
    plt.ylabel('f(x)',fontsize=14)
    plt.title('n = '+str(idx-1),fontsize=16)
    plt.grid()
    plt.tight_layout()
    
    plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap2_dichotomie_animation_'+str(idx)+'.png')
    
    if (f(a_n)*f(x_n))<0:
        b_n = x_n
        
    if (f(x_n)*f(b_n))<0:
        a_n = x_n
        
    x_n = (a_n+b_n)/2
