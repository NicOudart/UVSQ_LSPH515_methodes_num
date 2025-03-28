import numpy as np
import matplotlib.pyplot as plt

# def g(x):
    
#     return 2/x

# def g(x):
    
#     return 2*x-2/x

def g(x):
    
    return x/2+1/x

x = np.linspace(1.3,2.1,1000)

x_n = 2

#Loop:----------------------------------

for idx in range(0,4):
    
    plt.figure(idx)
    plt.plot(x,g(x),'r-')
    plt.plot(x,x,'k-')
    plt.plot([x_n,x_n],[x_n,g(x_n)],'k--')
    plt.plot([x_n,g(x_n)],[g(x_n),g(x_n)],'g--')
    plt.scatter(x_n,x_n,color='k',marker='o')
    plt.scatter(g(x_n),g(x_n),color='g',marker='o')
    plt.scatter(2**0.5,2**0.5,color='r',marker='o')
    plt.annotate('x_'+str(idx),(x_n,x_n),textcoords="offset points",xytext=(0,25),ha='center',fontsize=14)
    plt.annotate('x_'+str(idx+1),(g(x_n),g(x_n)),textcoords="offset points",xytext=(0,-30),ha='center',fontsize=14,color='g')
    plt.xlim([1.3,2.1])
    plt.ylim([1,2.5])
    plt.xlabel('x',fontsize=14)
    plt.ylabel('g(x)',fontsize=14)
    plt.title('n = '+str(idx),fontsize=16)
    plt.grid()
    plt.tight_layout()
    
    x_n = g(x_n)
    
    plt.savefig('C:/Users/oudart/Documents/Enseignements/Blog_PH515/Chap2_point_fixe_g3_animation_'+str(idx)+'.png')