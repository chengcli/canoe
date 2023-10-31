#! /usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
from netCDF4 import Dataset
from pylab import *

case = 'straka-test-main'
data = Dataset('%s.nc' % case, 'r')
time = data['time'][:]
x1 = data['x1'][:]/1.E3
x2 = data['x2'][:]/1.E3

tslice = [0., 300., 600., 900.]
Rd = 287.
cp = 1004.
p0 = 1.E5
Ts = 300.

# plot potential temperature
fig, axs = subplots(4, 1, figsize = (10, 12), sharex = True, sharey = True)
subplots_adjust(hspace = 0.08)
X, Y = meshgrid(x2, x1)
clines = hstack((arange(-17., 0.), arange(1., 5.)))

for i, t in enumerate(tslice):
  j = where(time >= t)[0][0]
  theta = data['theta'][j,:,:,0]
  theta -= Ts
  print('time = %.2f, theta min = %.2f, theta max = %.2f' % (time[j], theta.min(), theta.max()))
  ax = axs[i]
  ax.contour(X, Y, theta, clines, colors = 'k', linestyles = '-')
  #ax.contour(X, Y, theta.T, clines, colors = 'k')
  ax.set_xlim([0., 19.2])
  ax.set_ylim([0., 4.8])
  ax.text(14., 4.4, 'time = %.1f ' % time[j], fontsize = 12)
  ax.text(14., 4.0, '%.2f m' % (1.E3*(max(x2) - min(x2))/len(x2),),
    fontsize = 12)
  ax.text(14., 3.6, 'min T = %.2f K' % theta.min(), fontsize = 12)
  ax.text(14., 3.2, 'max T = %.2f K' % theta.max(), fontsize = 12)
  ax.set_ylabel('Z (km)', fontsize = 15)
ax.set_xlabel('X (km)', fontsize = 15)
savefig('2d.straka-theta.png' % case, bbox_inches = 'tight')
close()
