#! /usr/bin/env python3
from pylab import *
from pyathena import FileMode, io_wrapper, parameter_input
from pyharp import radiation

file = io_wrapper()
inp = parameter_input()

file.open("giants_ir.inp", FileMode.read)
inp.load_from_file(file)
file.close()

rad = radiation()
rad.load_all_radiation_bands(inp)
print(rad.get_num_bands())
for i in range(rad.get_num_bands()):
    band = rad.get_band(i)
    print(band.get_name())
