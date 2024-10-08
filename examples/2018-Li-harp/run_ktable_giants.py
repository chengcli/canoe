#! /usr/bin/env python3
import sys

sys.path.append("../python")
sys.path.append(".")

import canoe
from canoe import load_configure, find_resource
from canoe.harp import radiation_band
from atm_profile_utils import read_atm_profile
from multiprocessing import Pool, cpu_count
from typing import Tuple
from rfmlib import *


def run_ktable_one_band(bname: str):
    band = radiation_band(bname, config)

    nspec = band.get_num_specgrids()
    wmin, wmax = band.get_range()
    wres = (wmax - wmin) / (nspec - 1)

    # atmospheric properties
    num_layers = len(atm["HGT"])
    wav_grid = (wmin, wmax, wres)

    species = list(map(str, config[bname]["opacity"]))
    driver = create_rfm_driver(wav_grid, tem_grid, species, hitran_file)

    # write rfm atmosphere file to file
    write_rfm_atm(atm, rundir=bname)
    write_rfm_drv(driver, rundir=bname)

    # run rfm and write kcoeff file
    pwd = os.getcwd()
    run_rfm(rundir=bname)
    write_ktable(
        opacity_output + "-" + bname, species, atm, wav_grid, tem_grid, basedir=pwd
    )
    os.chdir(pwd)


if __name__ == "__main__":
    canoe.start()

    atm_profile = "new_jup_lat30_modelc.atm"
    opacity_config = "giants-op.yaml"
    opacity_output = "giants-lbl"
    tem_grid = (5, -50, 50)
    max_threads = cpu_count()

    hitran_file = find_resource("HITRAN2020.par")
    atm = read_atm_profile(atm_profile)
    config = load_configure(opacity_config)

    bnames = list(map(str, config["bands"]))
    for bname in bnames:
        run_ktable_one_band(bname)

    # pool = Pool(max_threads)
    # pool.map(run_ktable_one_band, bnames)
