#! /usr/bin/env python3
import sys, os, subprocess

sys.path.append("/Users/chengcli/Development/canoe/build/python")

hitran_file = "/Users/chengcli/Development/canoe/data/HITRAN2020.par"

from pyharp import radiation_band
from utilities import load_file
from collections import OrderedDict
from numpy import *


def check_file_exist(filename: str) -> bool:
    if not os.path.isfile(filename):
        print("File %s does not exist!" % filename)
        sys.exit(1)
    return True


# create rfm atmosphere file
def create_rfm_atm(absorbers: list, atm: dict) -> None:
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
def create_rfm_drv(driver: dict) -> None:
    print("# Creating rfm.drv ...")
    with open("rfm.drv", "w") as file:
        for sec in driver:
            if driver[sec] != None:
                file.write(sec + "\n")
                file.write(" " * 4 + driver[sec] + "\n")
    print("# rfm.drv written.")


def run_rfm() -> None:
    process = subprocess.Popen(
        ["./rfm.release"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )

    for line in iter(process.stdout.readline, b""):
        # decode the byte string and end='' to avoid double newlines
        print(line.decode(), end="")

    process.communicate()


def create_netcdf_input(
    bane_name: str,
    absorbers: list,
    wmin: float,
    wmax: float,
    wres: float,
    tnum: int,
    tmin: float,
    tmax: float,
    atm: dict,
) -> str:
    print(f"# Creating kcoeff.inp-{band_name} ...")
    fname = "kcoeff.inp-%s" % band_name
    with open("kcoeff.inp-%s" % band_name, "w") as file:
        file.write("# Molecular absorber\n")
        file.write("%d\n" % len(absorbers))
        file.write(" ".join(absorbers) + "\n")
        file.write("# Molecule data files\n")
        for ab in absorbers:
            file.write("%-40s\n" % ("./tab_" + ab.lower() + ".txt",))
        file.write("# Wavenumber range\n")
        file.write(
            "%-14.6g%-14.6g%-14.6g\n" % (wmin, wmax, int((wmax - wmin) / wres) + 1)
        )
        file.write("# Relative temperature range\n")
        file.write("%-14.6g%-14.6g%-14.6g\n" % (tmin, tmax, tnum))
        file.write("# Number of vertical levels\n")
        file.write("%d\n" % len(atm["TEM"]))
        file.write("# Temperature\n")
        for i in range(len(atm["TEM"])):
            file.write("%-14.6g" % atm["TEM"][-(i + 1)])
            if (i + 1) % 10 == 0:
                file.write("\n")
        if (i + 1) % 10 != 0:
            file.write("\n")
        file.write("# Pressure\n")
        for i in range(len(atm["PRE"])):
            file.write("%-14.6g" % atm["PRE"][-(i + 1)])
            if (i + 1) % 10 == 0:
                file.write("\n")
        if (i + 1) % 10 != 0:
            file.write("\n")
    print(f"# kcoeff.inp-{band_name} written.")
    return fname


def run_kcoeff(inpfile: str, ncfile: str) -> None:
    process = subprocess.Popen(
        ["./kcoeff.release", "-i", inpfile, "-o", ncfile],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    for line in iter(process.stdout.readline, b""):
        # decode the byte string and end='' to avoid double newlines
        print(line.decode(), end="")

    process.communicate()


if __name__ == "__main__":
    info = load_file("example_earth.yaml")
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
    atm = {}
    ps_mbar = 1.0e3
    Ts_K = 300.0
    grav_ms2 = 9.8
    mu_kgmol = 29.0e-3
    Rgas_JKmol = 8.314
    Rgas_Jkg = Rgas_JKmol / mu_kgmol
    Hscale_m = Rgas_Jkg * Ts_K / grav_ms2
    num_layers = 100

    # km
    atm["HGT"] = linspace(0, 20, num_layers, endpoint=False)
    atm["PRE"] = ps_mbar * exp(-atm["HGT"] * 1e3 / Hscale_m)
    atm["TEM"] = Ts_K * ones(num_layers)
    # ppmv
    atm["CO2"] = 400.0 * ones(num_layers)

    create_rfm_atm(absorbers, atm)

    # create temperature pertubation
    tem_pertub = (5, -20, 20)
    wave_grid = (wmin, wmax, wres)

    # create rfm driver file
    driver = OrderedDict(
        [
            ("*HDR", "Header for rfm"),
            ("*FLG", "TAB CTM"),
            ("*SPC", "%.4f %.4f %.4f" % wave_grid),
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
    run_rfm()
    inp = create_netcdf_input(band_name, absorbers, *wave_grid, *tem_pertub, atm)
    run_kcoeff(inp, f"kcoeff-{band_name}.nc")
