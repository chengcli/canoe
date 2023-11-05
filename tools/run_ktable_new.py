#! /usr/bin/env python3
import sys, os, subprocess

sys.path.append("/Users/chengcli/Development/canoe/build/python")

hitran_file = "/Users/chengcli/Development/canoe/data/HITRAN2020.par"

from pyharp import radiation_band
from utilities import load_file
from collections import OrderedDict
from numpy import *


def check_file_exist(filename) -> bool:
    if not os.path.isfile(filename):
        print("File %s does not exist!" % filename)
        sys.exit(1)
    return True


# create rfm atmosphere file
def create_rfm_atm(absorbers, atm) -> None:
    print("# Creating rfm.atm ...")
    num_layers = atm["HGT"].shape[0]
    num_absorbers = len(absorbers)
    with open("rfm.atm", "w") as file:
        file.write("%d\n" % num_layers)
        file.write("*HGT [km]\n")
        for i in range(num_layers):
            file.write("%.8g " % atm["HGT"][i])
        file.write("\n*PRE [mb]\n")
        for i in range(num_layers):
            file.write("%.8g " % atm["PRE"][i])
        file.write("\n*TEM [K]\n")
        for i in range(num_layers):
            file.write("%.8g " % atm["TEM"][i])
        for i in range(num_absorbers):
            name = absorbers[i]
            file.write("\n*" + name + " [ppmv]\n")
            for j in range(num_layers):
                file.write("%.8g " % atm[name][j])
        file.write("\n*END")
    print("# rfm.atm written.")


# create rfm driver file
def create_rfm_drv(driver) -> None:
    print("# Creating rfm.drv ...")
    with open("rfm.drv", "w") as file:
        for sec in driver:
            if driver[sec] != None:
                file.write(sec + "\n")
                print(driver[sec])
                file.write(" " * 4 + driver[sec] + "\n")
    print("# rfm.drv written.")


def run_rfm() -> None:
    print("# Running rfm ...")
    script = ["rfm", "rfm.drv"]
    out, err = subprocess.Popen(
        script, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    ).communicate()
    os.system("rfm")
    print("# rfm finished.")


if __name__ == "__main__":
    info = load_file("example_earth.yaml")
    opacity = info["opacity-sources"]
    band_name = str(info["bands"][0])

    band = radiation_band(band_name, info)

    nspec = band.get_num_spec_grids()
    wmin, wmax = band.get_range()
    wres = (wmax - wmin) / (nspec - 1)
    print(wmin, wmax, wres)

    num_absorbers = band.get_num_absorbers()
    absorbers = []
    for i in range(num_absorbers):
        absorbers.append(band.get_absorber(i).get_name())

    # create atmosphere dictionary
    atm = {}
    ps_pa = 1.0e5
    Ts_K = 300.0
    grav_ms2 = 9.8
    mu_kgmol = 29.0e-3
    Rgas_JKmol = 8.314
    Rgas_Jkg = Rgas_JKmol / mu_kgmol
    Hscale_m = Rgas_Jkg * Ts_K / grav_ms2
    num_layers = 100

    atm["HGT"] = linspace(0, 20, num_layers, endpoint=False) * 1.0e3
    atm["PRE"] = ps_pa * exp(-atm["HGT"] / Hscale_m)
    atm["TEM"] = Ts_K * ones(num_layers)
    atm["CO2"] = 400.0 * ones(num_layers)

    create_rfm_atm(absorbers, atm)

    # create temperature pertubation
    tem_pertub = (5, -20, 20)

    # create rfm driver file
    driver = OrderedDict(
        [
            ("*HDR", "Header for rfm"),
            ("*FLG", "TAB CTM"),
            ("*SPC", "%.4f %.4f %.4f" % (wmin, wmax, wres)),
            ("*GAS", " ".join(absorbers)),
            ("*ATM", "rfm.atm"),
            ("*DIM", "PLV \n    %d %.4f %.4f" % tem_pertub),
            ("*TAB", "tab_*.txt"),
            ("*HIT", hitran_file),
            ("*END", ""),
        ]
    )
    check_file_exist(hitran_file)
    create_rfm_drv(driver)
