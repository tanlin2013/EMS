import numpy as np
from scipy import integrate
from scipy.misc import derivative
import matplotlib.pyplot as plt
from scipy import inf
from scipy.optimize import fsolve
from sympy.solvers import solve
from sympy import Symbol
import sympy as sp
import time
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt

start = time.time()
class BH:
    def __init__(self,a,b,c,mu,dz):
        self.a=a
        self.b=b
        self.c=c
        self.mu=mu
        self.dz=dz        
        #self.dz0=dz0
    def partial_derivative(self,func,axis=0,point=[],n=1):
        args=point[:]
        def func_along(x):
            args[axis]=x
            return func(*args)
        return derivative(func_along,point[axis],self.dz,n)    
    
    def Ae(self,z):
        Ae=-self.a*np.log(self.b*z**2+1)
        return Ae 
    
    def phi(self,z):
        phi_p=lambda y: np.sqrt(-6*(derivative(self.Ae,y,self.dz,n=2)-derivative(self.Ae,y,self.dz)**2+2*derivative(self.Ae,y,self.dz)/y))
        phi=integrate.quad(phi_p,0.,z)
        return phi[0]
        
    def As(self,z):
        return self.Ae(z)+np.sqrt(1.0/6.0)*self.phi(z)
    
    def g(self,z,zH):
        def I(self,yi,yf):
            f=lambda y: y**3*np.exp(-3*self.Ae(y))
            I=integrate.quad(f,yi,yf)
            return I[0]
            
        def I2(self,yi,yf):
            f2=lambda y: y**3*np.exp(-3*self.Ae(y))*np.exp(self.c*y**2)
            I2=integrate.quad(f2,yi,yf)
            return I2[0] 
            
        def det(z,zH):
            vec=np.array([[I(self,0.,zH) , I2(self,0.,zH)],
                          [I(self,zH,z) , I2(self,zH,z)]])
            det=np.linalg.det(vec)
            return det
            
        g=1+1./I(self,0.,zH) \
        *(2*self.c*self.mu**2/((1-np.exp(self.c*zH**2))**2)*det(z,zH)-I(self,0.,z))
        return g
        
    def g2(self,z,zH):
        def I(self,yi,yf):
            f=lambda y: y**3*np.exp(-3*self.Ae(y))
            I=integrate.quad(f,yi,yf)
            return I[0]
            
        def I2(self,yi,yf):
            f=lambda y: y*np.exp(self.c*y**2)
            I2=integrate.quad(f,yi,yf)
            return I2[0]
            
        def I3(self,yi,yf):
            f=lambda x,y: y**3*np.exp(-3*self.Ae(y))*x*np.exp(self.c*x**2) 
            I3=integrate.dblquad(f,yi,yf,lambda x: 0.,lambda x: x)            
            return I3[0]
        
        g2p=lambda y: -y**3*np.exp(-3*self.Ae(y))/I(self,0.,zH) \
        *(1-(self.mu/I2(self,0.,zH))**2*(I2(self,0.,y)*I(self,0.,zH)-I3(self,0.,zH)))
        g2=integrate.quad(g2p,0.,z)
        return g2[0]+1.0
        
    def V(self,z,zH):
        V=-3*z**2*self.g(z)*np.exp(-2*self.Ae(z)) \
        *(derivative(self.Ae,z,self.dz,n=2) + 3*derivative(self.Ae,z,self.dz)**2 \
        +(1.5*self.partial_derivative(self.g,axis=0,point=[z,zH])/self.g(z,zH)-6./z)*derivative(self.Ae,z,self.dz) \
        -(1.5*self.partial_derivative(self.g,axis=0,point=[z,zH])/self.g(z,zH)-4./z)/z \
        +self.partial_derivative(self.g,axis=0,point=[z,zH],n=2)/(6*self.g(z,zH)))
        return V
        
    def T(self,zH):
        T=-self.partial_derivative(self.g,axis=0,point=[zH,zH])/(4.*np.pi)
        return T
    
    def sBH(self,z):
        s=np.exp(3*self.Ae(z))/(4.*z**3)
        return s
        
    def F(self,zH):
        f=lambda y: self.sBH(y)*derivative(self.T,y,self.dz)
        F=integrate.quad(f,zH,2.)
        return F[0]
        
    def sigma(self,z,zH):
        sigma=np.exp(2*(self.Ae(z)+np.sqrt(1.0/6.0)*self.phi(z)))*np.sqrt(self.g(z,zH))/(z**2)        
        return sigma
        
    def r(self,z0,zH):
        f=lambda y: 2.0/np.sqrt(self.g(y,zH)*((self.sigma(y,zH)/self.sigma(z0,zH))**2-1.0))
        r=integrate.quad(f,0.,z0)        
        return r[0]
        
    def Vm(self,z0,zH):
        f= lambda y: 2.0*(self.sigma(y,zH)/np.sqrt(self.g(y,zH)*(1.0-(self.sigma(z0,zH)/self.sigma(y,zH))**2))-(1.0+2.0*np.sqrt(6.0*self.a*self.b)*y)/y**2)
        #e= lambda x:-2.0*(1.0+2.0*np.sqrt(self.a*self.b)*x)/x**2
        #-2.0/z0+4.0*np.sqrt(self.a*self.b)*np.log(z0)+
        Vm=integrate.quad(f,0.01,z0)
        return Vm[0]-2.0/z0+4*np.sqrt(self.a*self.b)*np.log(z0)
    def C(self,z0):
        return -2.0/z0+4*np.sqrt(self.a*self.b)*np.log(z0)
    def intVm(self,z,z0,zH):
        return 2.0*(self.sigma(z,zH)/np.sqrt(self.g(z,zH)*(1.0-(self.sigma(z0,zH)/self.sigma(z,zH))**2))-(1.0+2.0*np.sqrt(6.0*self.a*self.b)*z)/z**2)
    def sVm(self,z0,zH):
        N=1000
        yi=0.001
        yf=z0-0.00001
        h=np.float((yf-yi)/N)
        s=np.float((self.intVm(yi,z0,zH)+self.intVm(yf,z0,zH)))
        for i in range (1,N/2+1):
            s+=4.0*np.float(self.intVm(yi+(2*i-1)*h,z0,zH))
        for i in range (1,N/2):
            s+=2.0*np.float(self.intVm(yi+2*i*h,z0,zH))
        return (h/3.0)*s    
   
    def numVm(self,z0,zH):
        return self.sVm(z0,zH)+self.C(z0)
    def sigmas(self,zm,zH):
        return self.sigma(zm,zH)
    def zm(self,zH):
        f=lambda y: self.partial_derivative(self.sigma,axis=0,point=[y,zH])
        return fsolve(f,1.0)[0]
    
    def Sth(self,zH,zt):
        Sth=
        return self.
############################################################
        
if __name__=='__main__':

    a=4.046
    b=0.01613
    c=0.227
    mu=0.12    
    dz=0.01 # spacing for numerical derivative       
    sim=BH(a,b,c,mu,dz)
################################
    zH=1.0
    numVmlist1=[] 
    z_interval1=np.linspace(0.001,zH,100)
    
    for z0 in z_interval1:              
        numVm=sim.numVm(z0,zH)
        numVmlist1.append(numVm)
        
    zH=3.0
    z_interval2=np.linspace(0.001,zH,100)   
    numVmlist2=[]   
    for z0 in z_interval2:       
        numVm=sim.numVm(z0,zH)
        numVmlist2.append(numVm)    
    zH=10.0
    z_interval3=np.linspace(0.001,zH,100)
    numVmlist3=[]   
    for z0 in z_interval3:
        numVm=sim.numVm(z0,zH)
        numVmlist3.append(numVm) 
    plt.plot(z_interval1,numVmlist1,'b-')
    plt.plot(z_interval1,numVmlist2,'r-')
    plt.plot(z_interval1,numVmlist3,'c-')
    plt.legend()
    plt.savefig("fig_Vz.eps", format="eps") 
    
###################################################    
    zH=1.0
    numVmlist1=[] ;rlist1=[]
    z_interval1=np.linspace(0.001,0.768,100)
    
    x_extra1=np.linspace(1.083,16.0,100)
    y_extra1=np.linspace(0.817,0.817,100)
    
    for z0 in z_interval1:      
        r=sim.r(z0,zH)
        rlist1.append(r)
        numVm=sim.numVm(z0,zH)
        numVmlist1.append(numVm)
        
    zH=3.0
    z_interval2=np.linspace(0.001,2.243,100)
    
    x_extra2=np.linspace(5.36,16.0,100)
    y_extra2=np.linspace(10.68,10.68,100)
    
    numVmlist2=[] ;rlist2=[]    
    for z0 in z_interval2:      
        r=sim.r(z0,zH)
        rlist2.append(r)
        numVm=sim.numVm(z0,zH)
        numVmlist2.append(numVm)    
    zH=10.0
    z_interval3=np.linspace(0.001,4.041,40)
    numVmlist3=[] ;rlist3=[]    
    for z0 in z_interval3:      
        r=sim.r(z0,zH)
        rlist3.append(r)
        numVm=sim.numVm(z0,zH)
        numVmlist3.append(numVm)
    
    plt.plot(rlist1,numVmlist1,'b-')
    plt.plot(rlist2,numVmlist2,'r-')   #r-z figure   
    plt.plot(rlist3,numVmlist3,'c-')
    plt.plot(x_extra1,y_extra1,'b-')
    plt.plot(x_extra2,y_extra2,'r-')
    
    plt.xlim(0,15)
    plt.ylim(-10,25)
    plt.ylabel(r'$V(GeV)$',family='Times New Roman',style='italic',rotation=0,labelpad=25) 
    plt.xlabel(r'$r(GeV^{-1})$',family='Times New Roman',style='italic')
    plt.text(9, -3, r'$T>T_{\mu}$', fontsize=15)    
    plt.text(9, 7, r'$T=T_{\mu}$', fontsize=15)
    plt.text(9, 13.3, r'$T<T_{\mu}$', fontsize=15)  
        
    
    plt.legend()
    plt.savefig("fig_Vr.eps", format="eps")    
##########################################################  
    zH=3.5
    z0_interval=np.linspace(0.0,zH,100)
    z_interval=np.linspace(0.0,1.0,100)
    r1list=[];r2list=[] 
    for z0 in z0_interval:
       r=sim.r(z0,zH)
       r1list.append(r)
       
    x_extra1=np.linspace(0,1,100)   
    y_extra1=np.linspace(7.3,7.3,100)
   
    zH=9.0
    z0_interval=np.linspace(0.0,3.909,25)
    for z0 in z0_interval:
       r=sim.r(z0,zH)
       r2list.append(r)       
    x_extra2=np.linspace(0.434,0.434,100)
    y_extra2=np.linspace(0,60.0,100)
    z_interval2=np.linspace(0.0,0.434,25)
    
    x_extra3=np.linspace(0.747,0.747,100)
    y_extra3=np.linspace(0,7.3,100)

    plt.plot(z_interval,r1list,'b-')
    plt.plot(z_interval2,r2list,'r-')
    plt.plot(x_extra2,y_extra2,'r--')
    plt.plot(x_extra1,y_extra1,'b--')
    plt.plot(x_extra3,y_extra3,'b--')
    
    plt.xlim(0,1)
    plt.ylim(0,35)
    
    plt.ylabel(r'$r(GeV^{-1})$',family='Times New Roman',style='italic',rotation=0,labelpad=25) 
    plt.xlabel(r'$z_0/z_H$',family='Times New Roman',style='italic')
    plt.text(-0.04, 7.2, r'$r_M$', fontsize=8)    
    plt.text(0.713, -1.5, r'$z_M/z_H$', fontsize=8)
    plt.text(0.4, -1.5, r'$z_m/z_H$', fontsize=8)
    plt.xticks([0,0.33,0.66,0.9,1.0])     
    plt.legend()
    plt.savefig("fig_rz.eps", format="eps")
#############################################################    
   
    
end = time.time()
elapsed = end - start
print "Time taken: ", elapsed, "seconds."
     
