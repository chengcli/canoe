#! /usr/bin/env python3
from pylab import *
from pyathena import FileMode, io_wrapper, parameter_input
from pycanoe import index_map
import pyharp
from pyharp import radiation

file = io_wrapper()
pin = parameter_input()

file.open("giants_ir.inp", FileMode.read)
pin.load_from_file(file)
file.close()

idx = index_map.init_from_athena_input(pin)

pyharp.init_index_map(pin)

rad = radiation()
rad.load_all_radiation_bands(pin)
print(rad.get_num_bands())
for i in range(rad.get_num_bands()):
    #band = rad.get_band(i)
    print(rad.get_band(i).get_name())
