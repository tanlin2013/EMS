import bg
import numpy as np
from scipy import integrate
from scipy.misc import derivative

def _partial_derivative(func,axis=0,point=[],n=1,dx=1e-2):
    args = point[:]
    def func_along(x):
        args[axis] = x
        return func(*args)
    return derivative(func_along,point[axis],dx,n)

def s_BH(model,zh):
    s_BH = bg.we(model.Ae,zh)**model.d / 4.
    return s_BH

def T_BH(model,zh,mu,dz=1e-2):
    T_BH = - _partial_derivative(bg.g,axis=1,point=[model,zh,zh,mu],dx=dz)/(4*np.pi)
    return T_BH

def F_BH(model,zh,mu,dz=1e-2,ulim=10):
    func = lambda y: s_BH(model,y)*_partial_derivative(T_BH,axis=1,point=[model,y,mu],dx=dz)
    F_BH, err = integrate.quad(func,zh,ulim)
    return F_BH

def sigma(model,z,zh,mu):
    sigma = bg.ws(model,z)**2 * np.sqrt(bg.g(model,z,zh,mu))
    return sigma
    
def Vqq(model,zh,mu,zt,eps=1e-2):
    func = lambda y: 2/(4*np.pi) * bg.ws(model,y)**2 * (1 - (sigma(model,zt,zh,mu)/sigma(model,y,zh,mu))**2)**(-0.5)
    Vqq, err = integrate.quad(func,eps,zt)
    return Vqq

def rqq(model,zh,mu,zt):
    func = lambda y: 2*(bg.g(model,z,zh,mu)*((sigma(model,y,zh,mu)/sigma(model,zt,zh,mu))**2-1))**(-0.5)
    rqq,err = integrate.quad(func,0,zt)
    return rqq
    
def S_TH(model,zh,mu,zt,dz=1e-2):
    S_TH = -_partial_derivative(Vqq,axis=1,point=[model,zh,mu,zt],dx=dz)/_partial_derivative(T_BH,axis=1,point=[model,z,zh,mu],dx=dz)
    #func = lambda y: ws(model,zt)**4/ws(model,z)**2 (1 - (sigma(model,zt,zh,mu,tol)/sigma(model,y,zh,mu,tol))**2)**(-1.5) 
           #*_partial_derivative(bg.g(model,zt,zh,mu,tol)/bg.g(model,z,zh,mu,tol),axis=2,point=[model,y,mu,tol],dx=tol)
    #S_TH, err = integrate.quad(func,0,zt) /_partial_derivative(bg.g,axis=1,point=[model,zh,zh,mu,dz],dx=dz)
    return S_TH

def F_EF(model,zh,ztmtol=1e-2):
    F_EF = T_BH()*
    return F_EF
