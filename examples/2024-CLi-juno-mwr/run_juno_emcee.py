#! /usr/bin/env python3
import numpy as np
import emcee
import sys, os
import matplotlib.pyplot as plt
import h5py

sys.path.append("../python")
sys.path.append(".")

from canoe import def_species, load_configure
from canoe.snap import def_thermo
from canoe.athena import Mesh, ParameterInput, Outputs, MeshBlock
# from canoe.harp import radiation_band, radiation

# forward operater
def run_RT_modify_atmos(mb: MeshBlock, 
                 adlnTdlnP: float=0., 
                 pmin: float=0. ,
                 pmax: float =0.
                 )-> np.array:
    # adlnTdlnP=0.0 ## set as insensitive
    mb.modify_dlnTdlnP(adlnTdlnP, pmin, pmax)
    # adlnNH3dlnP = 0#.25
    # mb.modify_dlnNH3dlnP(adlnNH3dlnP, pmin, pmax)

    # for k in range(mb.k_st, mb.k_ed + 1):
    #    for j in range(mb.j_st, mb.j_ed + 1):
    rad=mb.get_rad()
    rad.cal_radiance(mb, mb.k_st, mb.j_st)

    nb=rad.get_num_bands()
    tb=np.array([0.]*4*nb)

    for ib in range(nb):
        print(rad.get_band(ib))
        toa=rad.get_band(ib).get_toa()[0]
        tb[ib*4:ib*4+4]=toa
    return tb

def set_atmos_run_RT(qNH3: float, # ppmv 
                     T0: float=180.  # Kelvin
                     ):
    
    mb.construct_atmosphere(pin,qNH3,T0)    
    rad=mb.get_rad()
    rad.cal_radiance(mb, mb.k_st, mb.j_st)

    nb=rad.get_num_bands()
    tb=np.array([0.]*4*nb)

    for ib in range(nb):
        # print(rad.get_band(ib))
        toa=rad.get_band(ib).get_toa()[0]
        tb[ib*4:ib*4+4]=toa
    # print(tb[4:])
    return tb[4:]

# Define likelihood function
def ln_likelihood(theta, observations, observation_errors):
    nh3, temperature = theta

    simulations = set_atmos_run_RT(nh3,temperature)  # Use your forward operator here
    residuals = observations - simulations
    print(simulations)
    print(observations)
    # print(residuals)
    chi_squared = np.sum((residuals / observation_errors) ** 2)
    # print(chi_squared)
    return -0.5 * chi_squared

# Define priors for NH3 and temperature
def ln_prior(theta):
    nh3, temperature = theta

    nh3_mean = 300  # Mean value for NH3
    nh3_stddev = 100  # Standard deviation for NH3
    temperature_mean = 169  # Mean value for temperature
    temperature_stddev = 10  # Standard deviation for temperature   0.5%
    
    ln_prior_nh3 = -0.5 * ((nh3 - nh3_mean) / nh3_stddev)**2 - np.log(nh3_stddev * np.sqrt(2 * np.pi))
    ln_prior_temperature = -0.5 * ((temperature - temperature_mean) / temperature_stddev)**2 - np.log(temperature_stddev * np.sqrt(2 * np.pi))
    
    if (0 < nh3 < 1000):# and (100 < temperature < 200):
        return ln_prior_nh3 + ln_prior_temperature
    return -np.inf  # return negative infinity if parameters are outside allowed range

# Combine likelihood and prior to get posterior
def ln_posterior(theta, observations, observation_errors):
    prior = ln_prior(theta)
    if not np.isfinite(prior):
        return -np.inf
    return prior + ln_likelihood(theta, observations, observation_errors)


## main
if __name__ == "__main__":

    ##  extract TB observations from ZZ fitting results
    observations=np.zeros((20,))  
    obs=np.zeros((24,))  
    pj=51
    mu=np.cos(np.array([0.,15.,30.,45.])/180.*np.pi)
    print(mu)
    for ch in range(6):
        tb_file=h5py.File(f"/nfs/nuke/chengcli/JUNOMWR/zzhang/PJ{pj:02d}_Freq{ch}.h5","r")
        if ch==0:
            c0=tb_file['ModelTypeupdate1_MultiPJ_Mode1/Iter1/c0'][-1] ## the north polar 
            c1=tb_file['ModelTypeupdate1_MultiPJ_Mode1/Iter1/c1'][-1] ## the north polar 
            c2=tb_file['ModelTypeupdate1_MultiPJ_Mode1/Iter1/c2'][-1] ## the north polar 
        else:
            c0=tb_file['ModelTypeupdate1_MultiPJ_Mode3/Iter2/c0'][-1] ## the north polar 
            c1=tb_file['ModelTypeupdate1_MultiPJ_Mode3/Iter2/c1'][-1] ## the north polar 
            c2=tb_file['ModelTypeupdate1_MultiPJ_Mode3/Iter2/c2'][-1] ## the north polar 
        tb_file.close()

        # Xr=1.0 ## \mu >0.6
        obs[(ch)*4:(ch+1)*4]=c0-c1*5.*(1-mu)+c2/0.04*0.5*(mu-0.8)*(1-mu)
    ## discard CH1
    observations=obs[4:]
    print(observations)

    # [740.51932939 732.39178625 708.02076917 667.58562359 474.58510281
    # 469.42513666 454.00808555 428.59559118 338.13016122 335.65949356
    # 328.01197674 314.60534003 251.9730167  250.46642377 245.71888005
    # 237.15115289 194.47971955 193.67185714 191.10407859 186.40702401
    # 141.18445694 141.06723252 140.59821156 139.46693178]

    ## initialize Canoe
    global pin
    pin = ParameterInput()
    pin.load_from_file("juno_mwr.inp")

    vapors = pin.get_string("species", "vapor").split(", ")
    clouds = pin.get_string("species", "cloud").split(", ")
    tracers = pin.get_string("species", "tracer").split(", ")

    def_species(vapors=vapors, clouds=clouds, tracers=tracers)
    def_thermo(pin)
    
    config = load_configure("juno_mwr.yaml")
    # print(pin.get_real("problem", "qH2O.ppmv"))

    mesh = Mesh(pin)
    mesh.initialize(pin)

    global mb
    mb = mesh.meshblock(0)

##  run MCMC

    # Generate synthetic observations and errors (replace with your real data)
    # observations = np.random.normal(size=20)  # Replace with your observations
    observation_errors_stddev = 0.03 * observations   ## error = 5% 

    # Initialize walkers
    n_walkers = 5
    n_dimensions = 2  # nh3, temperature
    initial_guess = [200.0, 150.0]  # Initial guess for NH3 and temperature
    initial_guesses = [initial_guess + 10*np.random.randn(n_dimensions) for _ in range(n_walkers)]

    filename = "run_mcmc_5000.h5"
    backend = emcee.backends.HDFBackend(filename)
    backend.reset(n_walkers, n_dimensions)

    # Set up the sampler
    sampler = emcee.EnsembleSampler(n_walkers, n_dimensions, ln_posterior, args=(observations, observation_errors_stddev), backend=backend)

    # Run MCMC
    n_steps = 5000
    sampler.run_mcmc(initial_guesses, n_steps, progress=True)

    # Extract samples
    samples = sampler.get_chain()

    # Compute mean and standard deviation of samples
    mean_nh3 = np.mean(samples[:,:,0])
    std_nh3 = np.std(samples[:,:,0])
    mean_temperature = np.mean(samples[:,:,1])
    std_temperature = np.std(samples[:,:,1])

    print("NH3: Mean =", mean_nh3, "Standard Deviation =", std_nh3)
    print("Temperature: Mean =", mean_temperature, "Standard Deviation =", std_temperature)


    # Plot convergence diagnostics
    fig, axes = plt.subplots(2, figsize=(10, 7), sharex=True)
    for iwalk in range(n_walkers):
        axes[0].plot(range(n_steps),sampler.get_chain()[:,iwalk,0].T, alpha=0.4)
    axes[0].set_ylabel("NH3")
    for iwalk in range(n_walkers):
        axes[1].plot(range(n_steps),sampler.get_chain()[:,iwalk,1].T, alpha=0.4)
    axes[1].plot(sampler.get_chain()[:,:,1].T, alpha=0.4)
    axes[1].set_ylabel("Temperature")
    axes[1].set_xlabel("Step")
    plt.savefig("MCMC-5000Steps.png")

    # Flatten samples for posterior distribution plotting
    flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)

    # Plot posterior distributions
    labels = ["NH3", "Temperature"]
    fig, axes = plt.subplots(2, figsize=(10, 7), sharex=True)
    for i in range(2):
        ax = axes[i]
        ax.hist(flat_samples[:, i], bins=50, color="skyblue", histtype="stepfilled", alpha=0.7)
        ax.set_title(labels[i])
        ax.grid(True)
    # plt.tight_layout()
    plt.savefig("posterior-distributions-5000.png")
