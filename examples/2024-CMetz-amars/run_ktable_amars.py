#! /usr/bin/env python3
import sys, os

sys.path.append("../python")
sys.path.append(".")

from canoe import def_species, load_configure, find_resource
from canoe.harp import radiation_band
from netCDF4 import Dataset
from numpy import *
from rfmlib import *


# create reference atmosphere dictionary
# Edit this function to create your own atmosphere
def create_rfm_atmosphere(nlyr: int) -> dict:
    atm = {}

    data = Dataset("amars-main.nc", "r")
    # m
    atm["HGT"] = (data["x1f"][:-1] + data["x1f"][1:]) / 2.0

    # pa
    atm["PRE"] = data["press"][0, :, 0, 0]

    # K
    atm["TEM"] = data["temp"][0, :, 0, 0]

    # mole fraction
    atm["H2O"] = data["vapor1"][0, :, 0, 0] / 18.0 * 2.2

    # mole fraction
    atm["H2S"] = data["vapor2"][0, :, 0, 0] / 17.0 * 2.2

    # mole fraction
    atm["SO2"] = data["vapor3"][0, :, 0, 0] / 17.0 * 2.2

    # mole fraction
    atm["CO2"] = data["vapor4"][0, :, 0, 0] / 17.0 * 2.2

    return atm


if __name__ == "__main__":
    hitran_file = find_resource("HITRAN2020.par")
    def_species(vapors=["H2O", "H2S", "SO2"])

    config = load_configure("amars.yaml")

    for i in range(8):
        band = radiation_band(str(config["bands"][i]), config)

        nspec = band.get_num_specgrids()
        wmin, wmax = band.get_range()
        wres = (wmax - wmin) / (nspec - 1)

        absorbers = ["H2O", "H2S", "SO2", "CO2"]

        # atmospheric properties
        num_layers = 128
        wav_grid = (wmin, wmax, wres)
        tem_grid = (5, -50, 50)

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
        write_ktable(fname + "_B" + str(i + 1), absorbers, atm, wav_grid, tem_grid)
