# import emcee

# %config InlineBackend.figure_format = "retina"
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# reader = emcee.backends.HDFBackend("run_mcmc.h5")

# tau = reader.get_autocorr_time()
# burnin = int(2 * np.max(tau))
# thin = int(0.5 * np.min(tau))
# samples = reader.get_chain(discard=burnin, flat=True, thin=thin)
# log_prob_samples = reader.get_log_prob(discard=burnin, flat=True, thin=thin)

import h5py

h5 = h5py.File("run_mcmc_5000.h5", "r")
chain = h5["mcmc"]["chain"][:]
h5.close()

flattened_chain = chain.reshape(-1, 2)

# Create labels for the parameters
labels = ["qNH3 [ppmv]", "Temperature [K]"]

# Create the corner plot
fig, ax = plt.subplots(2, 1, figsize=(17, 10))
for iw in range(5):
    ax[0].plot(range(5000), chain[:, iw, 0])

ax[0].set_ylabel("qNH3 [ppmv]")
ax[0].set_xlim([0, 5000])

for iw in range(5):
    ax[1].plot(range(5000), chain[:, iw, 1])
ax[1].set_ylabel("Temperature [K]")
ax[1].set_xlim([0, 5000])
ax[1].set_xlabel("step")
plt.tight_layout()
# Show the plot
plt.savefig("emcee_inspect_step_5000.png")
