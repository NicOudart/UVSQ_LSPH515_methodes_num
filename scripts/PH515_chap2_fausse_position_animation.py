import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return (x**2)-2

x = np.linspace(-0.5,2.5,1000)

a_n = 0
b_n = 2

x_n = a_n

#Loop:----------------------------------

for idx in range(1,6):

    plt.figure(idx)
    plt.plot(x,f(x),'r-')
    plt.scatter(2**0.5,0,color='r',marker='o')
    
    q_n = (f(b_n)-f(a_n))/(b_n-a_n)
    
    plt.plot(x,q_n*x+(f(a_n)-q_n*a_n),'k-')
    plt.scatter([a_n,b_n],[f(a_n),f(b_n)],color='k',marker='o')
    plt.annotate('x_'+str(idx),(a_n,f(a_n)),textcoords="offset points",xytext=(0,25),ha='center',fontsize=14)
    plt.annotate('x_'+str(idx-1),(b_n,f(b_n)),textcoords="offset points",xytext=(0,-25),ha='center',fontsize=14)
    
    x_n = x_n-f(x_n)/q_n
    
    if (f(a_n)*f(x_n))<0:
        b_n = x_n
        
    if (f(x_n)*f(b_n))<0:
        a_n = x_n

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
    
    plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap2_fausse_position_animation_'+str(idx)+'.png')
    
    
