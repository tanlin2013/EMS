import EMS
import numpy as np
import matplotlib.pyplot as plt

#class heavy_quark: 

class light_quark:
    def __init__(self):
        self.a = 4.046
        self.b = 0.01613
        self.c = 0.227
        self.d = 3

    def Ae(self,y):
        Ae = -self.a*np.log(self.b*y**2+1)
        return Ae
    
    def f(self,y):
        f = np.exp(-self.c*y**2-self.Ae(y))
        return f
    
if __name__=='__main__':    
    
    model = light_quark()
    zs=np.linspace(0.,30.,80)
    zh=10.0
    mu=1.0
    
    gs = []
    for z in zs:
        g = EMS.bg.blackening_factor(model.d,z,zh,mu,model.f,model.Ae)
        gs.append(g)
        
    plt.figure()
    plt.plot(zs,gs,marker='o',linestyle='-')
    plt.xlim(0,12)
    plt.ylim(-0.5,1.2)
    plt.grid()
    plt.show()
