#! /usr/bin/env python3
import sys, os

sys.path.append("../python")
sys.path.append(".")

from pyharp import radiation_band, subscribe_species
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

    # ppmv
    atm["H2O"] = data["vapor1"][0, :, 0, 0] / 18.0 * 2.2 * 1.0e6

    # ppmv
    atm["NH3"] = data["vapor2"][0, :, 0, 0] / 17.0 * 2.2 * 1.0e6

    return atm


if __name__ == "__main__":
    subscribe_species(
        {
            "vapor": ["H2O", "NH3"],
        }
    )

    config = load_configure("jupiter_rt.yaml")
    band = radiation_band("B1", config, load_opacity=True)

    aH2O = band.get_absorber_by_name("H2O")
    aNH3 = band.get_absorber_by_name("NH3")
    aH2H2CIA = band.get_absorber_by_name("H2-H2-CIA")
    aH2HeCIA = band.get_absorber_by_name("H2-He-CIA")

    temp = 300.0
    xH2O = 0.003
    xNH3 = 0.0003
    v1 = 0.0
    v2 = 0.0
    v3 = 0.0
    pres = 10.0e5

    atm = [temp, xH2O, xNH3, v1, v2, v3, pres]
    print(aH2O.get_attenuation(1000.0, 1000.0, atm))
    print(aNH3.get_attenuation(1000.0, 1000.0, atm))
    print(aH2H2CIA.get_attenuation(1000.0, 1000.0, atm))
    print(aH2HeCIA.get_attenuation(1000.0, 1000.0, atm))

    num_layers = 100
    band.resize(num_layers)

    # atm = create_atmosphere(num_layers)
    # band.set_spectral_properties(atm)

    tau = band.get_tau()
    print(tau.shape)
