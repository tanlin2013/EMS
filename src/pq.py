import bg
import numpy as np
from scipy.misc import derivative

def _partial_derivative(func,axis=0,point=[],n=1,dx=1e-8):
    args = point[:]
    def func_along(x):
        args[axis] = x
        return func(*args)
    return derivative(func_along,point[axis],dx,n)

def s_BH(model,zh,dz=1e-8):
    s_BH = bg.we(model.Ae,zh)**model.d / 4.
    return s_BH

def T_BH(model,zh,mu,dz=1e-8):
    T_BH = - _partial_derivative(bg.blackening_factor,axis=1,point=[model,zh,zh,mu,dz],dx=dz)/(4*np.pi)
    return T_BH

def F_BH(model,zh,mu,dzh=1e-8):
    func = lambda y: s_BH(model,y)*_partial_derivative(T_BH,axis=1,point=[model,y,mu,dzh],dx=dz)
    F_BH, err = integrate.quad(func,zh,np.inf)
    return F_BH

def sigma(model,z,zh,mu,tol=1e-8):
    sigma = bg.ws(model,z)**2 * np.sqrt(bg.blackening_factor(model,z,zh,mu,tol))
    return sigma
    
def Vqq(model,zh,zt,mu):
    func = lambda y: 2/(4*np.pi) * bg.ws(model,y)**2 * (
        1 - (sigma(model,zt,zh,mu)/sigma(model,y,zh,mu))**2)**(-0.5)
    Vqq, err = integrate.quad(func,0,zt)
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
