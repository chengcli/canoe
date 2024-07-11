#! /usr/bin/env python3
import sys, os

sys.path.append("../python")
sys.path.append(".")

from pyharp import radiation_band, subscribe_species
from utilities import load_configure, find_resource
from netCDF4 import Dataset
from numpy import *
from rfmlib import *


# create reference atmosphere dictionary
# Edit this function to create your own atmosphere
def create_rfm_atmosphere(nlyr: int) -> dict:
    atm = {}

    data = Dataset("hjupiter-main.nc", "r")
    # km
    atm["HGT"] = (data["x1f"][:-1] + data["x1f"][1:]) / (2.0 * 1.0e3)

    # mbar
    atm["PRE"] = data["press"][0, :, 0, 0] / 100.0

    # K
    atm["TEM"] = data["temp"][0, :, 0, 0]

    # ppmv
    atm["H2"] = 1.0e6

    return atm


if __name__ == "__main__":
    config = load_configure("hjupiter.yaml")
    band = radiation_band("B1", config)

    nspec = band.get_num_spec_grids()
    wmin, wmax = band.get_range()
    wres = (wmax - wmin) / (nspec - 1)

    absorbers = ["H2O", "NH3"]

    # atmospheric properties
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

    fname = band.get_absorber_by_name("H2O").get_opacity_file()
    base_name = os.path.basename(fname)
    fname, _ = os.path.splitext(base_name)
    write_ktable(fname, absorbers, atm, wav_grid, tem_grid)
