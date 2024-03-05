from constants import *
import numpy as np
from numpy import sin, cos, cosh, sqrt, pi, arctan, tanh, sinh
from gw_functions import phiv, omg, get_M
from hypmik3pn import get_u_hat
from scipy.integrate import cumtrapz
import antenna_pattern as ap

    

def get_hyp_waveform(M,q,et,n0,t,inc,distance,phi0='None'):
    if phi0=='None':
        phi0=0


    η=q/(1+q)**2
    Time=M*tsun
    dis=M*dsun
    scale=distance/dis
    
    
    l=n0*t
    x0=(M*tsun*n0)**(2/3)
    u=get_u_hat(l,et)

    
    phi=phiv(η,et,u,x0,order=3)
    r1=(-1+et*cosh(u))/x0;z=1/r1
    z=1/r1
    rt=(et*sinh(u)/(-1+et*cosh(u)))*sqrt(x0)
    phit=((et**2-1)**(1/2)/(-1+et*cosh(u))**2)*x0**(3/2)
    X=r1*cos(phi)
    Y=r1*sin(phi)
    
    hp_arr=(-η*(sin(inc)**2*(z-r1**2*phit**2-rt**2)+(1+cos(inc)**2)*((z
        +r1**2*phit**2-rt**2)*cos(2*phi)+2*r1*rt*phit*sin(2*phi))))
    hx_arr=(-2*η*cos(inc)*((z+r1**2*phit**2-rt**2)*sin(2*phi)-2*r1*rt*phit*cos(2*phi)))
    Hp=hp_arr/scale; Hx=hx_arr/scale
    
    


    return Hp,Hx





def cal_sp_sx_A(M,q,et,n0,t,D_GW,inc):
    η=q/(1+q)**2
    Time=M*tsun
    dis=M*dsun
    scale=D_GW/dis
    l=n0*t
    u=get_u_hat(l,et)
    Ci=cos(inc)
    Si=sin(inc)
    x0=(M*tsun*n0)**(2/3)
    ω=omg(η,et,u,x0)
    dd=(et*cosh(u)-1)
    PH=sqrt(et**2-1)*(cosh(2*u)-et*cosh(u))/dd
    QH=(et+(et**2-2)*cosh(u))*sinh(u)/dd
    RH=et*sinh(u)
    sp=η*x0/n0*((Ci**2+1)*(PH*sin(2*ω)-QH*cos(2*ω))+Si**2*RH)
    sx=η*x0/n0*2*Ci*(-PH*cos(2*ω)-QH*sin(2*ω))
    return sp/scale, sx/scale


def cal_sp_sx(M,q,et,n,tarr,D_GW,inc):
    
    η=1/(1+q)**2
    h=get_hyp_waveform(M, q, et, n,tarr, inc, D_GW)
    s_arr = np.array([cumtrapz(h[i], x = tarr, initial=0) for i in range(len(h))])
    
    
    return s_arr




from enterprise.signals.signal_base import function as enterprise_function, PTA
from scipy.interpolate import CubicSpline

@enterprise_function
def hyp_pta_res(toas,
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
    cons_terms=None,
    interp_steps=1000
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
    b           is the impact parameter of the GW source in total mass
    e0          is the eccentricity of the GW source
    log10_z     is the log10 cosmological redshift of the GW source
    tref        is the fiducial time in s in SSB frame
    interp_steps is the number of samples used for interpolating the PTA signal
    """
    order = 3
    M = 10**log10_M # Solar mass
    S = 10**log10_S
    n= 10**log10_n
    dis=M*dsun
    Time=M*tsun
    
    

    ts = toas - tref

    # ti, tf, tzs in seconds, in source frame
    ti = min(ts)
    tf = max(ts)
    Tspan=tf-ti
    eta=q/(1+q)**2
    
    x0=(M*n*tsun)**(2/3)
    
    D_GW = dis*x0*eta/(S*n)
    
    tz_arr = np.linspace(ti, tf, interp_steps)
    delta_t_arr = (tz_arr[1]-tz_arr[0])

    
    
   

    inc = np.arccos(cos_inc)

    gwra = gwphi
    gwdec = np.arcsin(cos_gwtheta)

    psrra = phi
    psrdec = np.pi/2 - theta
    
    # if method=='num':
    #     s_arr=cal_sp_sx(M,q,e0,b,tz_arr,D_GW,inc)
    # elif method=='analytic':
    #     sA=cal_sp_sx_A(M,q,e0,b,tz_arr,D_GW,inc)
    #     spA=sA[0]-sA[0][0];sxA=sA[1]-sA[1][0]
    #     s_arr=[spA,sxA]
    sA=cal_sp_sx_A(M,q,e0,n,tz_arr,D_GW,inc)
    
    if(cons_terms is None):
        
        hp0,hx0=get_hyp_waveform(M,q,e0,n,ts[0],inc,D_GW)
        dt=tz_arr-ti;sp0,sx0=(hp0*dt,hx0*dt)
    
    elif(cons_terms is True):
        (sp0,sx0)=(0,0)
    
    spA=sA[0]-sA[0][0]-sp0;sxA=sA[1]-sA[1][0]-sx0
    s_arr_=[spA,sxA]
    
    cosmu, Fp, Fx = ap.antenna_pattern(gwra, gwdec, psrra, psrdec)
    
    c2psi = np.cos(2*psi)
    s2psi = np.sin(2*psi)
    Rpsi = np.array([[c2psi, -s2psi],
                     [s2psi, c2psi]])
    
    s_arr = np.dot([Fp,Fx], np.dot(Rpsi, s_arr_))
    
    s_spline = CubicSpline(tz_arr, s_arr)
    s = s_spline(ts)

    return s
    
    
    
    
    
    
