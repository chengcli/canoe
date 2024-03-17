# Purpose: Library for RFM calculations

import os, subprocess
from collections import OrderedDict


# create rfm driver file
def create_rfm_driver(
    wav_grid: tuple[float, float, float],
    tem_grid: tuple[float, float, float],
    absorbers: list[str],
    hitran_file: str,
) -> dict:
    driver = OrderedDict(
        [
            ("*HDR", "Header for rfm"),
            ("*FLG", "TAB CTM"),
            ("*SPC", "%.4f %.4f %.4f" % wav_grid),
            ("*GAS", " ".join(absorbers)),
            ("*ATM", "rfm.atm"),
            ("*DIM", "PLV \n    %d %.4f %.4f" % tem_grid),
            ("*TAB", "tab_*.txt"),
            ("*HIT", hitran_file),
            ("*END", ""),
        ]
    )
    return driver


# write rfm atmosphere file
def write_rfm_atm(atm: dict) -> None:
    print("# Creating rfm.atm ...")
    num_layers = atm["HGT"].shape[0]
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
        for name, val in atm.items():
            if name in ["HGT", "PRE", "TEM"]:
                continue
            file.write("\n*" + name + " [ppmv]\n")
            for j in range(num_layers):
                file.write("%.8g " % val[j])
        file.write("\n*END")
    print("# rfm.atm written.")


# write rfm driver file
def write_rfm_drv(driver: dict) -> None:
    print("# Creating rfm.drv ...")
    with open("rfm.drv", "w") as file:
        for sec in driver:
            if driver[sec] != None:
                file.write(sec + "\n")
                file.write(" " * 4 + driver[sec] + "\n")
    print("# rfm.drv written.")


# run rfm
def run_rfm() -> None:
    process = subprocess.Popen(
        ["./rfm.release"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )

    for line in iter(process.stdout.readline, b""):
        # decode the byte string and end='' to avoid double newlines
        print(line.decode(), end="")

    process.communicate()


# create netcdf input file
def create_netcdf_input(
    fname: str,
    absorbers: list,
    atm: dict,
    wmin: float,
    wmax: float,
    wres: float,
    tnum: int,
    tmin: float,
    tmax: float,
) -> str:
    print(f"# Creating {fname}.inp ...")
    with open(f"{fname}.inp", "w") as file:
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
    print(f"# {fname}.inp written.")
    return fname


# write kcoeff table to netcdf file
def write_ktable(
    fname: str,
    absorbers: list[str],
    atm: dict,
    wav_grid: tuple[float, float, float],
    tem_grid: tuple[float, float, float],
) -> None:
    inpfile = create_netcdf_input(fname, absorbers, atm, *wav_grid, *tem_grid)

    process = subprocess.Popen(
        ["./kcoeff.release", "-i", inpfile, "-o", f"{fname}.nc"],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    for line in iter(process.stdout.readline, b""):
        # decode the byte string and end='' to avoid double newlines
        print(line.decode(), end="")

    process.communicate()
    print(f"# {fname}.nc written.")
