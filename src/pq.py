from scipy.misc import derivative

def _partial_derivative(func,axis=0,point=[],n=1,dx=1e-8):
    args = point[:]
    def func_along(x):
        args[axis] = x
        return func(*args)
    return derivative(func_along,point[axis],dx,n)
    
def T_BH(g,zh,dz=1e-8):
    T_BH = - derivative(g,zh,dz)/(4*np.pi)
    return T_BH
    
#def sigma():
#    sigma = 
#    return sigma
