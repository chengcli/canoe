#! /usr/bin/env python3
from netCDF4 import Dataset
from numpy import sum, sqrt
import sys

case1 = sys.argv[1]
case2 = sys.argv[2]

data1 = Dataset(case1, "r")
var1 = data1["vort"][-1, :, :, :]

data2 = Dataset(case2, "r")
var2 = data2["vort"][-1, :, :, :]

diff = sqrt(sum((var2 - var1) * (var2 - var1)))

if diff < 1.0:
    print("### polar vortex case passed. ###")
else:
    raise ValueError("ERROR: polar vortex failed. L2-norm is %.2g" % diff)
