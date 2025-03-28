import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return (x**2)-2

x = np.linspace(-0.5,2.5,1000)

x_n = 0
x_n_old = 2

#Loop:----------------------------------

for idx in range(1,6):

    plt.figure(idx)
    plt.plot(x,f(x),'r-')
    plt.scatter(2**0.5,0,color='r',marker='o')
    
    q_n = (f(x_n)-f(x_n_old))/(x_n-x_n_old)
    
    plt.plot(x,q_n*x+(f(x_n)-q_n*x_n),'k-')
    plt.scatter([x_n_old,x_n],[f(x_n_old),f(x_n)],color='k',marker='o')
    plt.annotate('x_'+str(idx),(x_n,f(x_n)),textcoords="offset points",xytext=(0,25),ha='center',fontsize=14)
    plt.annotate('x_'+str(idx-1),(x_n_old,f(x_n_old)),textcoords="offset points",xytext=(0,-25),ha='center',fontsize=14)
    
    x_n_old = x_n
    x_n = x_n - f(x_n)/q_n

    plt.plot([x_n,x_n],[0,f(x_n)],'g--')
    plt.scatter(x_n,f(x_n),color='g',marker='o')
    plt.annotate('x_'+str(idx+1),(x_n,f(x_n)),textcoords="offset points",xytext=(0,-25),ha='center',fontsize=14,color='g')
    plt.xlim([-0.5,2.5])
    plt.ylim([-3,3])
    plt.xlabel('x',fontsize=14)
    plt.ylabel('f(x)',fontsize=14)
    plt.title('n = '+str(idx),fontsize=16)
    plt.grid()
    plt.tight_layout()
    
    plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap2_secante_animation_'+str(idx)+'.png')
    
    
