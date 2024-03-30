#! /usr/bin/env python3
import sys, os

sys.path.append("../python")
sys.path.append(".")

from pyharp import radiation_band, init_thermo
from utilities import load_configure
from numpy import linspace, ones, exp
from netCDF4 import Dataset
from pylab import *


# create atmosphere dictionary
def create_atmosphere(nlyr: int) -> dict:
    atm = {}

    data = Dataset("jupiter1d-main.nc", "r")
    # km
    atm["HGT"] = (data["x1f"][:-1] + data["x1f"][1:]) / (2.0 * 1.0e3)

    # mbar
    atm["PRE"] = data["press"][0, :, 0, 0] / 100.0

    # K
    atm["TEM"] = data["temp"][0, :, 0, 0]

    return atm


if __name__ == "__main__":
    config = load_configure("hjupiter.yaml")
    init_thermo(config["thermodynamics"])
    band = radiation_band("ck", config, load_opacity=True)

    ab = band.get_absorber_by_name("premix")
    print(ab)

    temp = 300.0
    pres = 10.0e5
    v1 = 0.0
    v2 = 0.0
    v3 = 0.0

    air = [temp, v1, v2, v3, pres]
    print(ab.get_attenuation(1000.0, 1000.0, air))
