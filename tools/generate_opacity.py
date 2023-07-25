#! /usr/bin/env python3
from pylab import *
from pyathena import FileMode, io_wrapper, parameter_input
from pyharp import radiation, init_index_map

file = io_wrapper()
pin = parameter_input()

file.open("giants_ir.inp", FileMode.read)
pin.load_from_file(file)
file.close()

init_index_map(pin)

rad = radiation()
rad.load_all_radiation_bands(pin)
for i in range(rad.get_num_bands()):
    band = rad.get_band(i)
    print(band.get_name())
