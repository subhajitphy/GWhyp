from mpmath import besselk
import numpy as np
import scipy
from scipy import optimize
import matplotlib.pyplot as plt
from getx import get_x


def Fomg(eta,b1,e,omg):
    if type(omg) is np.ndarray:
        return np.array([Fomg(eta,b1,e,omg1) for omg1 in omg])
    b=b1
    a=b/np.sqrt(e**2-1)
    u=omg*e*a**(3/2)
    p=1j*u/e
    bk1 = besselk(p+1,u)
    bk0 = besselk(p,u)
    fac=32/5*(eta/a)**2*(p/u**2)**2*np.exp(-1j*np.pi*p)
    sed=(u**2*(p**2+u**2+1)*(p**2+u**2)*bk1**2
    -2*u*((p-3/2)*u**2+p*(p-1)**2)*(p**2+u**2)*bk0*bk1
    +2*(u**6/2+(2*p**2-3/2*p+1/6)*u**4+(5/2*p**4-7/2*p**3+p**2)*u**2+p**4*(p-1)**2)*bk0**2  )
    return float(np.abs(fac*sed))

def get_max(eta,b1,e1):
    x0=get_x(e1,eta,b1,0)[0]
    n=x0**(3/2)
    val=optimize.fmin(lambda omg:-Fomg(eta,b1,e1,omg),n,  xtol=1e-10, ftol=1e-10,disp= False)[0]
    return val
    #return float(np.abs(val))



