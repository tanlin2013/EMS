import numpy as np
from scipy import integrate

class integral_set:
    def __init__(self,d,f,Ae,tol):
        self.d = d
        self.f = f
        self.Ae = Ae
        self.tol = tol
        
    def I1(self,yi,yf):
        func = lambda y: we(y,self.Ae)**(2-self.d)/self.f(y) 
        I1, err = integrate.quad(func,yi,yf)
        return I1
        
    def I2(self,yi,yf):
        func = lambda y: we(y,self.Ae)**(-self.d)
        I2, err = integrate.quad(func,yi,yf)
        return I2
        
    def I12(self,x2i,x2f):
        func = lambda x1, x2: we(x2,self.Ae)**(-self.d) * we(x1,self.Ae)**(2-self.d)/self.f(x1)
        I12, err = integrate.dblquad(func, x2i, x2f, lambda x1: 0, lambda x1: x1)
        return I12
    
    def det(self,z,zh):
        mat = np.array([[self.I2(0,zh), self.I12(0,zh)],
                        [self.I2(0,z), self.I12(0,z)]])
        det = np.linalg.det(mat)
        return det
    
def blackening_factor(d,z,zh,mu,f,Ae,tol=1e-8):
    Iset = integral_set(d,f,Ae,tol)
    g = 1 - Iset.I2(0,z)/Iset.I2(0,zh) + mu**2 * (
            Iset.det(z,zh)/(Iset.I2(0,zh) * Iset.I1(0,zh)**2))
    return g

def we(z,Ae): 
    we = np.exp(Ae(z))/z
    return we
    
