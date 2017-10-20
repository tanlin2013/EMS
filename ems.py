import EMS
import numpy as np
import matplotlib.pyplot as plt

def Ae(y):
    Ae = -a*np.log(b*y**2+1)
    return Ae
    
def f(y):
    f = np.exp(-c*y**2-Ae(y))
    return f

if __name__=='__main__':
    
    a=4.046
    b=0.01613
    c=0.227
    d=3
    zs=np.linspace(0.,30.,80)
    zh=10.0
    mu=1.0
    
    gs = []
    for z in zs:
        g = EMS.bg.blackening_factor(d,z,zh,mu,f,Ae)
        gs.append(g)
        
    plt.figure()
    plt.plot(zs,gs,marker='o',linestyle='-')
    plt.xlim(0,12)
    plt.ylim(-0.5,1.2)
    plt.grid()
    plt.show()
