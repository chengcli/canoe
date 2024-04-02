# import emcee
import corner
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

h5=h5py.File('run_mcmc_5000.h5', 'r')
chain = h5['mcmc']['chain'][180:]
h5.close()

flattened_chain = chain.reshape(-1, 2)

# Create labels for the parameters
labels = ["qNH3 [ppmv]","Temperature [K]"]

# Create the corner plot
fig, ax = plt.subplots(2, 2, figsize=(10, 10))
corner.hist2d(flattened_chain[:, 0], flattened_chain[:, 1], ax=ax[1, 0])
ax[1, 0].set_xlabel("qNH3 [ppmv]")
ax[1, 0].set_ylabel("Temperature [K]")
ax[1, 0].set_xlim([0,700])

# Plot histograms for each parameter
x = np.linspace(1, 700, 700)
for i in range(2):
    ax[i, i].hist(flattened_chain[:, i], bins=30, color='blue', alpha=0.7, density=True)
    ax[i, i].set_xlabel(labels[i])
    ax[i, i].set_ylabel('PDF')
    means=np.mean(flattened_chain[:, i])
    stdev=np.std(flattened_chain[:, i])
    ax[i, i].axvline(means, color='b', linestyle='--', label=f'posterior: ({means:6.2f}, {stdev:5.2f}')

    # Plot prior distribution
    mean, stddev = [(300, 100), (169, 10)][i]

    prior = norm.pdf(x, mean, stddev)
    ax[i, i].plot(x, prior, color='red', linestyle='--', label=f'Prior: norm({mean}, {stddev})')
    ax[i, i].legend()
    x = np.linspace(131, 200, 70)

ax[0, 0].set_xlim([0,700])

ax[0, 1].axis('off')  # Disable the axis
ax[0, 1].set_visible(False)  # Hide the axis
# Show the plot
plt.tight_layout()
# Show the plot
plt.savefig("emcee_cornerplot_5000.png")