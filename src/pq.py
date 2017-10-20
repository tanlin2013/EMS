import bg
from scipy.misc import derivative

def _partial_derivative(func,axis=0,point=[],n=1,dx=1e-8):
    args = point[:]
    def func_along(x):
        args[axis] = x
        return func(*args)
    return derivative(func_along,point[axis],dx,n)

def s_BH(d,Ae,zh,dz=1e-8):
    s_BH = np.exp(3*bg.we(zh,Ae))/zh**3
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
    func = 
    Vqq, err = integrate.quad
    return Vqq

def rqq():
    func =
    rqq,err = integrate.quad
    return rqq
    
