import bg
import numpy as np
from scipy.misc import derivative

def _partial_derivative(func,axis=0,point=[],n=1,dx=1e-8):
    args = point[:]
    def func_along(x):
        args[axis] = x
        return func(*args)
    return derivative(func_along,point[axis],dx,n)

def s_BH(d,zh,Ae,dz=1e-8):
    s_BH = bg.we(zh,Ae)**d/4.
    return s_BH

def T_BH(g,zh,dz=1e-8):
    T_BH = - derivative(g,zh,dz)/(4*np.pi)
    return T_BH

def F_BH(g,zh,dz=1e-8):
    func = lambda y: s_BH()*
    F_BH, err = integrate.quad(func,yi,yf)
    return F_BH

def sigma(z,Ae):
    sigma = ws(z)**2*np.sqrt(g(z,Ae,dz))
    return sigma
    
def Vqq():
    func = 2/(4*np.pi)*(ws^2*(1-sigma(zt)^2)/sigma(z)^2)^(-1/2))
    Vqq, err = integrate.quad
    return Vqq

def rqq():
    func = 2*(g(z,zh,dz)*(sigma(z,Ae)^2)/sigma(zt)^2-1)^(-1/2)
    rqq,err = integrate.quad
    return rqq
    
def S_TH():
    func = 
    S_TH,err = integrate.quad
    return S_TH

def F_EF():
    F_EF = 
    return F_EF
