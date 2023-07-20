#! /usr/bin/env python3
from netCDF4 import Dataset
from numpy import sum, sqrt
import sys

case1 = sys.argv[1]
case2 = sys.argv[2]

data1 = Dataset(case1, 'r')
var1 = data1['rho'][-1,:,:,0]

data2 = Dataset(case2, 'r')
var2 = data2['rho'][-1,:,:,0]

diff = sqrt(sum((var2 - var1)*(var2 - var1)))

if diff < 1.E-6:
  print('### shallow water test passed. ###')
else:
  raise ValueError('ERROR: shallow water test failed. L2-norm is %.2g' % diff)
