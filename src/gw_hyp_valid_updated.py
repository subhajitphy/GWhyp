#Import libraries
from constants import *
import numpy as np
from numpy import sin, cos, cosh, sqrt, pi, sinh
from gw_functions import omg, rx
from hypmik3pn import get_u_hat
import antenna_pattern as ap
from enterprise.signals.signal_base import function as enterprise_function



#Analytic expressions for s+ and sx
def cal_sp_sx_A(M,q,et,n0,t,D_GW,inc,phi0):
    '''
    This function has the following changes:
    1. phi0 is now included as a parameter explicitly. It can take any value between 0 to pi.

    '''
    η=q/(1+q)**2
    Time=M*tsun
    dis=M*dsun
    scale=D_GW/dis
    l=n0*t
    u=get_u_hat(l,et)
    Ci=cos(inc)
    Si=sin(inc)
    C2p=cos(2*phi0)
    S2p=sin(2*phi0) 
    x0=(M*tsun*n0)**(2/3)
    ω=omg(η,et,u,x0)
    C2omg=cos(2*ω)
    S2omg=sin(2*ω)
    dd=(et*cosh(u)-1)
    PH=sqrt(et**2-1)*(cosh(2*u)-et*cosh(u))/dd
    QH=(et+(et**2-2)*cosh(u))*sinh(u)/dd
    RH=et*sinh(u)
    sp=η*x0/n0*((Ci**2+1)*S2p*(PH*C2omg+QH*S2omg) + (Ci**2+1)*C2p*(PH*S2omg-QH*C2omg) + Si**2*RH)
    sx=η*x0/n0*2*Ci*(-C2p*(PH*C2omg+QH*S2omg) + S2p*(PH*S2omg-QH*C2omg))
    return sp/scale, sx/scale



#Enterprise compatible function for residuals (phi0 is now a parameter and there is no interpolation anymore)
@enterprise_function
def hyp_pta_res_updated(toas,
    theta,
    phi,
    cos_gwtheta,
    gwphi,
    psi,
    cos_inc,
    log10_M,
    q,
    log10_n,
    e0,
    log10_S,
    tref,
    phi0,
    cons_off=None
):
    """
    Compute the PTA signal due to a hyperbolic encounter.
    toas        are pulsar toas in s in SSB frame
    theta       is pulsar zenith angle in rad
    phi         is pulsar RA in rad
    cos_gwtheta is cos zenith angle of the GW source
    gwphi       is the RA of the GW source in rad
    psi         is the GW polarization angle in rad
    cos_inc     is the cos inclination of the GW source
    log10_S     is the log10 amplitude of residuals of the GW source in sec
    q           is the mass ratio of the GW source
    log10_n     is the log10 of the characteristic frequency
    e0          is the eccentricity of the GW source
    tref        is the fiducial time in s in SSB frame
    phi0        is the orbital phase at t=tref
    cons_off    is whether there is a constant offset or not (None > Constant offset removed)
    ===========================================================================================

    This function has the following changes:
    1. phi0 is now included as a parameter explicitly. It can take any value between 0 to pi.
    2. Steps involving any sort of interpolation are removed since the model is completely analytic.
    3. The block involving removal of a linear-in-time trend is removed. The only option exists whether to remove a constant DC offset or not.

    """

    order = 3
    M = 10**log10_M # Solar mass
    S = 10**log10_S
    n= 10**log10_n
    dis=M*dsun
    Time=M*tsun
    
    
    ts = toas - tref

    eta=q/(1+q)**2   
    x0=(M*n*tsun)**(2/3)
    
    if rx(eta,e0,0,x0,3)<10:
        raise ValueError("rmin< 10M")
    
        
    D_GW = dis*x0*eta/(S*n)
    inc = np.arccos(cos_inc)

    gwra = gwphi
    gwdec = np.arcsin(cos_gwtheta)

    psrra = phi
    psrdec = np.pi/2 - theta
    
    sA=cal_sp_sx_A(M,q,e0,n,ts,D_GW,inc, phi0)
    
    if(cons_off is None):
        spA=sA[0]-sA[0][0];sxA=sA[1]-sA[1][0]
        s_arr_=[spA,sxA]

    elif(cons_off is True):
        spA=sA[0];sxA=sA[1]
        s_arr_=[spA,sxA]
    
    cosmu, Fp, Fx = ap.antenna_pattern(gwra, gwdec, psrra, psrdec)
    
    c2psi = np.cos(2*psi)
    s2psi = np.sin(2*psi)
    Rpsi = np.array([[c2psi, -s2psi],
                     [s2psi, c2psi]])
    
    s_arr = np.dot([Fp,Fx], np.dot(Rpsi, s_arr_))

    return s_arr