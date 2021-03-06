import pq
import numpy as np
from scipy import integrate
from scipy.misc import derivative

class _integral_set:
    def __init__(self,model,tol):
        self.d = model.d
        self.f = model.f
        self.Ae = model.Ae
        self.tol = tol

    def I1(self,yi,yf):
        func = lambda y: we(self.Ae,y)**(2-self.d)/self.f(y) 
        I1, err = integrate.quad(func,yi,yf)
        return I1
        
    def I2(self,yi,yf):
        func = lambda y: we(self.Ae,y)**(-self.d)
        I2, err = integrate.quad(func,yi,yf)
        return I2
        
    def I12(self,x2i,x2f):
        func = lambda x1, x2: we(self.Ae,x2)**(-self.d) * we(self.Ae,x1)**(2-self.d)/self.f(x1)
        I12, err = integrate.dblquad(func, x2i, x2f, lambda x1: 0, lambda x1: x1)
        return I12
    
    def det(self,z,zh):
        mat = np.array([[self.I2(0,zh), self.I12(0,zh)],
                        [self.I2(0,z), self.I12(0,z)]])
        det = np.linalg.det(mat)
        return det
    
def g(model,z,zh,mu,tol=1e-8):
    Iset = _integral_set(model,tol)
    g = 1 - Iset.I2(0,z)/Iset.I2(0,zh) + mu**2 * (
            Iset.det(z,zh)/(Iset.I2(0,zh) * Iset.I1(0,zh)**2))
    return g

def we(Ae,z): 
    we = np.exp(Ae(z))/z
    return we

def _dAe(Ae,z,dz=1e-2):
    dAe = derivative(Ae,z,dz)
    return dAe
    
def _d2Ae(Ae,z,dz=1e-2):
    d2Ae = derivative(Ae,z,dz,n=2)
    return d2Ae

def _dwe(Ae,z,dz=1e-2):
    dwe = pq._partial_derivative(we,axis=1,point=[Ae,z],dx=dz)
    return dwe
    
def _d2we(Ae,z,dz=1e-2):
    d2we = pq._partial_derivative(we,axis=1,point=[Ae,z],n=2,dx=dz)
    return d2we
    
def phi(model,z,eps=1e-4):
    func = lambda y: np.sqrt(-2*model.d * (_d2we(model.Ae,y)/we(model.Ae,y) - 2*(_dwe(model.Ae,y)/we(model.Ae,y))**2))
    #func = lambda y: np.sqrt(-2*model.d * (_d2Ae(model.Ae,y) - _dAe(model.Ae,y)**2 + 2*_dAe(model.Ae,y)/y))
    phi, err = integrate.quad(func,eps,z) 
    return phi

def As(model,z,eps=1e-4):
    As = model.Ae(z) + np.sqrt(1/6.)*phi(model,z,eps)
    return As

def ws(model,z,eps=1e-4):
    ws = np.exp(As(model,z,eps))/z
    return ws
