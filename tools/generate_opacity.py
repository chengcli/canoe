#! /usr/bin/env python3
import sys, os

sys.path.append("/Users/chengcli/Development/canoe/build/python")

from pyathena import FileMode, io_wrapper, parameter_input
from pyharp import radiation, init_index_map
from pycanoe import AirParcel

file = io_wrapper()
pin = parameter_input()

file.open("example_giants_ir.inp", FileMode.read)
pin.load_from_file(file)
file.close()

init_index_map(pin)
rad = radiation()
rad.load_all_radiation_bands(pin)

if air.hydro().flags["OWNDATA"]:
    print("The array owns its data.")
else:
    print("The array does not own its data.")
    # air.hydro[IPR] = 1.0

for i in range(rad.get_num_bands()):
    band = rad.get_band(i)
    name = band.get_name()
    wmin = band.get_wavenumber_min()
    wmax = band.get_wavenumber_max()
    wres = band.get_wavenumber_res()
    cia_list, hitran_list = [], []
    # print(band.get_num_absorbers())
    for j in range(band.get_num_absorbers()):
        if band.get_absorber(j).get_category() == "cia":
            cia_list.append(band.get_absorber(j).get_name())
        if band.get_absorber(j).get_category() == "hitran":
            hitran_list.append(band.get_absorber(j).get_name())
    print(name, wmin, wmax, wres, cia_list, hitran_list)
