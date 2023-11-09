#! /usr/bin/env python3
from canoe import *
from pyharp import radiation_band, subscribe_species
from utilities import load_file
from collections import OrderedDict
from numpy import *


# create atmosphere dictionary
def create_atmosphere(nlyr: int) -> dict:
    atm = {}
    ps_mbar = 1.0e3
    Ts_K = 300.0
    grav_ms2 = 9.8
    mu_kgmol = 29.0e-3
    Rgas_JKmol = 8.314
    Rgas_Jkg = Rgas_JKmol / mu_kgmol
    Hscale_km = Rgas_Jkg * Ts_K / grav_ms2 * 1.0e-3

    # km
    atm["HGT"] = linspace(0, 20, nlyr, endpoint=False)
    atm["PRE"] = ps_mbar * exp(-atm["HGT"] / Hscale_km)
    atm["TEM"] = Ts_K * ones(nlyr)
    # ppmv
    atm["CO2"] = 400.0 * ones(nlyr)
    # ppmv
    atm["H2O"] = 4000.0 * exp(-atm["HGT"] / (Hscale_km / 5.0))
    atm["O2"] = (1e6 - atm["CO2"] - atm["H2O"]) * 0.21
    atm["N2"] = (1e6 - atm["CO2"] - atm["H2O"]) * 0.78
    atm["Ar"] = (1e6 - atm["CO2"] - atm["H2O"]) * 0.01

    return atm


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

    for ab in band.absorbers():
        ab.load_opacity_from_file("kcoeff-B1.nc")

    num_layers = 100
    band.resize(num_layers)
    band.resize_solver(num_layers)

    atm = create_atmosphere(num_layers)
    band.set_spectral_properties(atm)

    toa = band.cal_radiance().get_toa()
    tau = band.get_tau()
    print(toa, tau.shape)
