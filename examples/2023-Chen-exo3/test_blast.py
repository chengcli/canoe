#! /usr/bin/env python3
from netCDF4 import Dataset
from numpy import sum, sqrt
import sys

case1 = sys.argv[1]
case2 = sys.argv[2]

data1 = Dataset(case1, "r")
rho1 = data1["rho"][-1, :, :, :]

data2 = Dataset(case2, "r")
rho2 = data2["rho"][-1, :, :, :]

diff = sqrt(sum((rho2 - rho1) * (rho2 - rho1)))

if diff < 10.0:
    print("### blast test passed. ###")
else:
    raise ValueError("ERROR: blas test failed. L2-norm is %.2g" % diff)
