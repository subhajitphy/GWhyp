from enterprise.signals.deterministic_signals import Deterministic
import sys
home="/mnt/hdd/Research/IPTA/"
sys.path.append(f'{home}/packages/GW_hyp/')
from gw_hyp import hyp_pta_res
from enterprise.signals.parameter import Uniform
import numpy as np

# For eDR3

tmin=4359621128.1630945
tmax=5153379820.877011

def gwhyp_1psr_block(
    tref=Uniform(tmin,tmax)("gwhyp_tref"),
    cos_gwtheta= Uniform(-1,1)("gwhyp_cgwt"),
    gwphi= Uniform(0,2*np.pi)("gwhyp_gwphi"),
    psi= Uniform(0,np.pi)("gwhyp_psi"),
    cos_inc= Uniform(-1,1)("gwhyp_cos_inc"),
    log10_M= Uniform(7,10)("gwhyp_log10_M"),
    q= Uniform(0.1,1)("gwhyp_q"),
    log10_n= Uniform(-10,-7.5)("gwhyp_log10_n"),
    e0= Uniform(1.08,2)("gwhyp_e0"),
    log10_S= Uniform(-10, -6)("gwhyp_log10_S"),
    name="gwhyp",
):
    return Deterministic(hyp_pta_res(
        cos_gwtheta=cos_gwtheta,
        gwphi=gwphi,
        psi=psi,
        cos_inc=cos_inc,
        log10_M=log10_M,
        q=q,
        log10_n=log10_n,
        e0=e0,
        log10_S=log10_S,
        tref=tref
             
    ),name=name)
