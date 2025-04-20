import numpy as np

N0 = 8000 #Paramètre de forme (m3/mm)
R = 5 #Taux de pluie (mm/h)

def f(D):
    
    return N0*np.exp(-(4.1*R**-0.21)*D)*D**6

D = np.linspace(0,6,1000)

#Methode de Gauss--------------------------------------------------------------

def gauss(f,a,b,y_i,v_i):
    
    #Initialiser la somme de la formule de quadrature :
    somme = 0
    
    #Boucles sur les points / poids :
    for i in range(len(y_i)):
        
        #Sommer l'évaluation de f au point choisi, pondérée par le poids choisi :
        somme += v_i[i]*f(y_i[i]*(b-a)/2+((a+b)/2))
    
    #Calculer l'estimation de l'aire sous la courbe :
    aire = somme*(b-a)/2
    
    return aire

print("Gauss d'ordre 0 = "+str(gauss(f,1,3,[0],[2])))
print("Gauss d'ordre 1 = "+str(gauss(f,1,3,[-(1/3)**0.5,(1/3)**0.5],[1,1])))
print("Gauss d'ordre 2 = "+str(gauss(f,1,3,[0,-(3/5)**0.5,(3/5)**0.5],[8/9,5/9,5/9])))
print("Gauss d'ordre 3 = "+str(gauss(f,1,3,[(3/7-2/7*(6/5)**0.5)**0.5,-(3/7-2/7*(6/5)**0.5)**0.5,(3/7+2/7*(6/5)**0.5)**0.5,-(3/7+2/7*(6/5)**0.5)**0.5],[(18+30**0.5)/36,(18+30**0.5)/36,(18-30**0.5)/36,(18-30**0.5)/36])))