#! /usr/bin/env python3
import matplotlib

matplotlib.use("Agg")
from netCDF4 import Dataset
from pylab import *

# case = 'straka-0320a'
case = "straka-b1-main"
data = Dataset("%s.nc" % case, "r")
time = data["time"][:]
x1 = data["x1"][:] / 1.0e3
x2 = data["x2"][:] / 1.0e3

Rd = 287.0
cp = 1004.0
p0 = 1.0e5
Ts = 300.0

# Potential temperature
# fig, axs = subplots(4, 1, figsize = (10, 12), sharex = True, sharey = True)
# subplots_adjust(hspace = 0.08)
X, Y = meshgrid(x2, x1)
clines = hstack((arange(-17.0, 0.0), arange(1.0, 5.0)))

for i in range(len(time)):
    figure(1, figsize=(12, 4))
    ax = axes()
    pres = data["press"][i, :, :, 0]
    rho = data["rho"][i, :, :, 0]
    temp = pres / (rho * Rd)
    pt = temp * (p0 / pres) ** (Rd / cp) - Ts
    print(data["press"].shape)
    print("time = %.2f, pt min = %.2f, pt max = %.2f" % (time[i], pt.min(), pt.max()))
    # ax.contourf(X, Y, pt, clines, cmap = 'Blues_r')
    ax.contour(X, Y, pt, clines, colors="k", linestyles="-")
    # ax.contour(X, Y, pt, clines, colors = 'k')
    ax.set_xlim([0.0, 19.2])
    ax.set_ylim([0.0, 4.8])
    ax.set_title("time = %.1f " % time[i], fontsize=12)
    # ax.text(14., 4.4, 'time = %.1f ' % time[i], fontsize = 12)
    # ax.text(14., 4.0, '%.2f m' % (1.E3*(max(x2) - min(x2))/len(x2),), fontsize = 12)
    # ax.text(14., 3.6, 'min T = %.2f K' % pt.min(), fontsize = 12)
    # ax.text(14., 3.2, 'max T = %.2f K' % pt.max(), fontsize = 12)
    savefig("%s-pt-%03d.png" % (case, i), bbox_inches="tight")
    close()
