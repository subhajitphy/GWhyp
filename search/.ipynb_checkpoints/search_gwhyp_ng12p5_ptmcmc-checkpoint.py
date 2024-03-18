import numpy as np
import enterprise.signals.parameter as parameter
import pickle
import json
import glob
from enterprise.signals import white_signals
import sys,os
from enterprise.signals import utils
from enterprise.pulsar import Pulsar
from enterprise.signals.parameter import Uniform
from enterprise.signals.gp_signals import MarginalizingTimingModel
from enterprise.signals.white_signals import MeasurementNoise
from enterprise.signals.signal_base import PTA
from enterprise.signals import selections
from enterprise.signals import gp_signals
from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc
from enterprise_extensions.sampler import JumpProposal as JP
from dynesty import plotting as dyplot
from matplotlib import pyplot as plt

import time

start=time.time()
####################################################################################################
#home="/mnt/hdd/Research/IPTA/"
#home="/data/home/prabu/InPTA/subhajit/"
home="/mnt/lustre/ibs/nrgw_school/astro/subhajit/work_dir/"
pgloc=f'{home}/packages/'
#pkls=glob.glob(f"{home}/pickles/*pkl")
datadir_=home+"datadir/NG12.5/"
datadir=datadir_
#datadir=datadir_+"GWHYP/Injected_NG/prpgwbs_Sm6nm8e1.16_ic/"
#######################################################################################################

#pkls=glob.glob(f"{datadir}/pickles/*pkl")
#psrs=[pickle.load(open(pkls[i], 'rb')) for i in range(len(pkls))]
psrlist=np.loadtxt(f"{datadir}/psrlist_44.txt",dtype="str")
parfiles=[sorted(glob.glob(f"{datadir}/{psr}*par"))[0] for psr in psrlist]
timfiles=[sorted(glob.glob(f"{datadir}/{psr}*tim"))[0] for psr in psrlist]

psrs = [Pulsar(par, tim, ephem="DE438") for par, tim in zip(parfiles, timfiles)]

psrlist = [psr.name for psr in psrs]

nfile=f'{datadir}/channelized_12p5yr_v3_full_noisedict.json'
with open(nfile, "r") as f:
    noisedict = json.load(f)
    
tmin = np.min([p.toas.min() for p in psrs])
tmax = np.max([p.toas.max() for p in psrs])
Tspan = tmax - tmin


name = "gwhyp"
priors={
 'cos_gwtheta': Uniform(-1,1)('cos_gwtheta'),
 'gwphi': Uniform(0,2*np.pi)('gwphi'),
 'psi': Uniform(0,np.pi)(f"{name}_psi"),
 'cos_inc': Uniform(-1,1)(f"{name}_cos_inc"),
 'log10_M': Uniform(7,10)(f"{name}_log10_M"),
 'q': Uniform(0.1,1)(f"{name}_q"),
 'log10_n': Uniform(-10,-6.5)(f"{name}_log10_n"),
 'e0': Uniform(1.08,2)(f"{name}_e0"),
 'log10_S': Uniform(-10, -5)(f"{name}_log10_S"),
 'tref': Uniform(tmin,tmax)("tref"),
 'cons_terms':True
}


from enterprise.signals.deterministic_signals import Deterministic
sys.path.append(f'{pgloc}/GW_hyp/')
from gw_hyp_valid import hyp_pta_res

wf=Deterministic(hyp_pta_res(**priors),name=name)



selection = selections.Selection(selections.by_backend)


# white noise parameters
efac = parameter.Constant()
equad = parameter.Constant()
ecorr = parameter.Constant()  # we'll set these later with the params dictionary

# red noise parameters
log10_A = parameter.Uniform(-20, -11)
gamma = parameter.Uniform(0, 7)

# log10_A = parameter.Constant()
# gamma = parameter.Constant()

# GW parameters (initialize with names here to use parameters in common across pulsars)
# white noise
ef = white_signals.MeasurementNoise(efac=efac, selection=selection)
eq = white_signals.TNEquadNoise(log10_tnequad=equad, selection=selection)
ec = white_signals.EcorrKernelNoise(log10_ecorr=ecorr, selection=selection)

# red noise (powerlaw with 30 frequencies)
pl = utils.powerlaw(log10_A=log10_A, gamma=gamma)
rn = gp_signals.FourierBasisGP(spectrum=pl, components=30, Tspan=Tspan)
tm = gp_signals.TimingModel(use_svd=True)

log10_A_gw = parameter.Uniform(-20, -11)("gwb_log10_A")
gamma_gw = parameter.Uniform(0, 7)("gwb_gamma")
    
    
cpl = utils.powerlaw(log10_A=log10_A_gw, gamma=gamma_gw)
gw = gp_signals.FourierBasisGP(spectrum=cpl, components=5, Tspan=Tspan, name="gwb")


ss=ef + eq + ec +  tm +  wf+gw+rn


models = [ss(p) for p in psrs]
pta =PTA(models)

pta.set_default_params(noisedict)



params = pta.param_names
ndim = len(params)

np.savetxt("pars.txt",pta.param_names, fmt="%s")
# groups = [list(np.arange(0, ndim))]
groups = []


# create parameter groups for the red noise parameters
rnpsrs = [p.split('_')[0] for p in params if 'red_noise_log10_A' in p]
for psr in rnpsrs:
    groups.extend([[params.index(psr + '_red_noise_gamma'), params.index(psr + '_red_noise_log10_A')]])    

if 'gwb_gamma' in params:
    groups.extend([[params.index('gwb_gamma'), params.index('gwb_log10_A')]])


x0=[]
for i in range(len(params)):
    if params[i] in noisedict:
        x0.append(noisedict[params[i]])
    else:
        x0.append(pta.params[i].sample())

x0 = np.hstack(x0)

x0 = np.hstack([p.sample() for p in pta.params])
ndim = len(x0)

cov = np.diag(np.ones(ndim) * 0.01**2)

chaindir='chains/opengwhyp/'

def gwhyp_target_likelihood_my(pta):
    def gwhyp_target_likelihood_fn(params):
        param_map = pta.map_params(params)
        try:
            lnlike = pta.get_lnlikelihood(param_map)
        except ValueError as err:
            print(err.args[0])
            lnlike = -np.inf
        return lnlike

    return gwhyp_target_likelihood_fn

get_lnlikelihood = gwhyp_target_likelihood_my(pta)

def likelihood(x):
    return get_lnlikelihood(x)

sampler = ptmcmc(ndim,
                 get_lnlikelihood,
                 pta.get_lnprior,
                 cov,
                 groups=groups,
                 resume=True,
                 outDir=chaindir)


emph_dir="/mnt/lustre/ibs/nrgw_school/astro/subhajit/work_dir/GWHYP.search/30Sep23/emph_dist/"

empirical_distr=pickle.load(open(f'{emph_dir}/empirical_distributions_5f.pkl', 'rb'))

jp = JP(pta, empirical_distr=empirical_distr)

sampler.addProposalToCycle(jp.draw_from_empirical_distr, 30)
sampler.addProposalToCycle(jp.draw_from_prior, 10)

# draw from ewf priors



# draw from gwb priors
gwb_params = [x for x in pta.param_names if "gwb" in x]
for param in gwb_params:
    sampler.addProposalToCycle(jp.draw_from_par_prior(param), 2)


Niter = int(5e6)
#x0 = np.hstack([p.sample() for p in pta.params])

sampler.sample(x0, Niter, SCAMweight=20, AMweight=25, DEweight=15)



end=time.time()

print("Execution time= %s min" % str((end-start)/60))


