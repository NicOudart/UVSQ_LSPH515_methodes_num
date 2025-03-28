import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return (x**2)-2

def f_derivee(x):

	return 2*x

x = np.linspace(1.2,2.75,1000)

x_n = 2.5

#Loop:----------------------------------

for idx in range(1,6):
    
    plt.figure(idx)
    plt.plot(x,f(x),'r-')
    plt.scatter(2**0.5,0,color='r',marker='o')
    
    q_n = f_derivee(x_n)
    
    plt.plot(x,q_n*x+(f(x_n)-q_n*x_n),'k-')
    plt.scatter(x_n,f(x_n),color='k',marker='o')
    plt.annotate('x_'+str(idx-1),(x_n,f(x_n)),textcoords="offset points",xytext=(0,-25),ha='center',fontsize=14)
    
    x_n = x_n - f(x_n)/q_n

    plt.plot([x_n,x_n],[0,f(x_n)],'g--')
    plt.scatter(x_n,f(x_n),color='g',marker='o')
    plt.annotate('x_'+str(idx),(x_n,f(x_n)),textcoords="offset points",xytext=(0,25),ha='center',fontsize=14,color='g')
    plt.xlim([1.2,2.75])
    plt.ylim([-0.5,4.5])
    plt.xlabel('x',fontsize=14)
    plt.ylabel('f(x)',fontsize=14)
    plt.title('n = '+str(idx-1),fontsize=16)
    plt.grid()
    plt.tight_layout()
    
    plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap2_newton_animation_'+str(idx-1)+'.png')
    
    
