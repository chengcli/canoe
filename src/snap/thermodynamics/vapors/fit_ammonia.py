#! /usr/bin/env python2.7
from scipy.stats import linregress
from pylab import *

data = genfromtxt('ammonia.dat')
temp = data[:,0] + 273.15
hfv = data[:,4]
a, b, r, p, e = linregress(temp, hfv)
print 'cpv = %.2f kJ/kg' % a

hfl = data[:,3]
a, b, r, p, e = linregress(temp, hfl)
print 'cpl = %.2f kJ/kg' % a

hfl = data[:,3]
a, b, r, p, e = linregress(temp, hfv - hfl)
print 'delta_cp = %.2f kJ/kg' % -a
print 'mu = %.2f kJ/kg' % b
