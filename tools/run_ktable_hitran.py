#! /usr/bin/env python3
import sys

sys.path.append("../python")

from pyharp import radiation_band, subscribe_species
from utilities import load_configure, find_resource
from numpy import *
from rfmlib import *


# create reference atmosphere dictionary
# Edit this function to create your own atmosphere
def create_ref_atmosphere(nlyr: int) -> dict:
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
    atm["NH3"] = 400.0 * ones(nlyr)
    # ppmv
    atm["H2O"] = 4000.0 * exp(-atm["HGT"] / (Hscale_km / 5.0))

    return atm


if __name__ == "__main__":
    # hitran_file = f"/home/chengcli/Model/canoe/data/HITRAN2020.par"
    hitran_file = find_resource("HITRAN2020.par")
    subscribe_species(
        {
            "vapor": ["H2O", "NH3"],
        }
    )

    config = load_configure("example_jupiter_rt.yaml")
    band = radiation_band("H2O-rotational", config)

    nspec = band.get_num_spec_grids()
    wmin, wmax = band.get_range()
    wres = (wmax - wmin) / (nspec - 1)

    absorbers = ["H2O", "NH3"]

    # atmospheric properties
    num_layers = 100
    wav_grid = (wmin, wmax, wres)
    tem_grid = (5, -20, 20)

    atm = create_ref_atmosphere(num_layers)
    driver = create_rfm_driver(wav_grid, tem_grid, absorbers, hitran_file)

    # write rfm atmosphere file to file
    write_rfm_atm(atm)
    write_rfm_drv(driver)

    # run rfm and write kcoeff file
    run_rfm()
    write_ktable(fname, absorbers, atm, wav_grid, tem_grid)
