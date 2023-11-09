#! /usr/bin/env python3
from canoe import *
from pyharp import radiation_band, subscribe_species
from utilities import load_file
from collections import OrderedDict
from run_ktable_simple import create_atmosphere
from numpy import *

if __name__ == "__main__":
    subscribe_species(
        {
            "vapor": ["H2O"],
            "tracer": ["CO2", "O3"],
        }
    )

    info = load_file("example_earth_rt.yaml")

    band_name = str(info["bands"][0])
    band = radiation_band(band_name, info)

    num_absorbers = band.get_num_absorbers()
    for i in range(num_absorbers):
        ab = band.get_absorber(i)
        ab.load_opacity_from_file("kcoeff-B1.nc")

    num_layers = 100
    band.resize(num_layers)
    band.resize_solver(num_layers)

    atm = create_atmosphere(num_layers)
    band.set_spectral_properties(atm)

    # print(band.btau)
    band.cal_band_radiance()
