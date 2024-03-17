#! /usr/bin/env python3
import sys, os

sys.path.append("../python")
sys.path.append(".")

from pyharp import radiation_band, subscribe_species
from utilities import load_configure
from numpy import linspace, ones, exp
from pylab import *


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
    atm["TEM"] = (
        Ts_K * ones(nlyr) * (atm["HGT"][-1] - atm["HGT"] + 10.0) / atm["HGT"][-1]
    )
    # atm["TEM"] = Ts_K * ones(nlyr)
    # ppmv
    atm["H2O"] = 4000.0 * exp(-atm["HGT"] / (Hscale_km / 5.0))

    return atm


if __name__ == "__main__":
    subscribe_species(
        {
            "vapor": ["H2O", "NH3"],
        }
    )

    config = load_configure("jupiter_rt.yaml")
    band = radiation_band("H2O-rotational", config, load_opacity=True)

    aH2O = band.get_absorber_by_name("H2O")
    aNH3 = band.get_absorber_by_name("NH3")
    aH2H2CIA = band.get_absorber_by_name("H2-H2-CIA")
    aH2HeCIA = band.get_absorber_by_name("H2-H2-CIA")

    exit()

    atm = [300.0, 0.01, 0.0, 0.0, 0.0, 1.0e5]
    k = ab.get_attenuation(130.0, 130.0, atm)
    print(k)

    num_layers = 100
    band.resize(num_layers)

    print("a1")

    atm = create_atmosphere(num_layers)
    band.set_spectral_properties(atm)

    print("a2")

    toa = band.cal_radiance().get_toa()
    tau = band.get_tau()

    print("a3")

    wavenumber = linspace(130, 260, int((260 - 130) / 0.01) + 1)

    print(wavenumber)
    print(len(wavenumber))
    print(toa.shape)
