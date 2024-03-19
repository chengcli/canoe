#! /usr/bin/env python3
import sys, os

sys.path.append("../python")
sys.path.append(".")

from pyharp import radiation_band
from utilities import load_configure, find_resource
from numpy import *
from rfmlib import *


# create atmosphere dictionary
def create_rfm_atmosphere(nlyr: int) -> dict:
    atm = {}
    Ps = 1.0e5
    Ts = 300.0
    grav = 9.8
    mu = 29.0e-3
    Rgas = 8.314 / mu
    Hscale = Rgas * Ts / grav

    # km
    atm["HGT"] = linspace(0, 20, nlyr, endpoint=False)
    atm["PRE"] = Ps * exp(-atm["HGT"] / Hscale)
    atm["TEM"] = Ts * ones(nlyr)
    # ppmv
    atm["CO2"] = 4.0e-4 * ones(nlyr)
    # ppmv
    atm["H2O"] = 4.0e-3 * exp(-atm["HGT"] / (Hscale / 5.0))
    atm["O2"] = (1.0 - atm["CO2"] - atm["H2O"]) * 0.21
    atm["N2"] = (1.0 - atm["CO2"] - atm["H2O"]) * 0.78
    atm["Ar"] = (1.0 - atm["CO2"] - atm["H2O"]) * 0.01

    return atm


if __name__ == "__main__":
    hitran_file = find_resource("HITRAN2020.par")

    info = load_configure("example_earth.yaml")
    opacity = info["opacity-sources"]
    band_name = str(info["bands"][0])

    band = radiation_band(band_name, info)

    nspec = band.get_num_spec_grids()
    wmin, wmax = band.get_range()
    wres = (wmax - wmin) / (nspec - 1)

    num_absorbers = band.get_num_absorbers()
    absorbers = []
    for i in range(num_absorbers):
        absorbers.append(band.get_absorber(i).get_name())

    # create atmosphere dictionary
    num_layers = 100
    wav_grid = (wmin, wmax, wres)
    tem_grid = (5, -20, 20)

    atm = create_rfm_atmosphere(num_layers)
    driver = create_rfm_driver(wav_grid, tem_grid, absorbers, hitran_file)

    # write rfm atmosphere file to file
    write_rfm_atm(atm)
    write_rfm_drv(driver)

    # run rfm and write kcoeff file
    run_rfm()
    write_ktable(f"kcoeff-{band_name}.nc", absorbers, atm, wav_grid, tem_grid)
