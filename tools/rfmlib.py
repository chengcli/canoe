# Purpose: Library for RFM calculations

import os, subprocess, shutil
import numpy as np
from collections import OrderedDict
from typing import List, Tuple, Dict

def create_rfm_driver(
    wav_grid: Tuple[float, float, float],
    tem_grid: Tuple[int, float, float],
    absorbers: List[str],
    hitran_file: str,
) -> Dict[str, str]:
    """
    Create a RFM driver file.

    Parameters
    ----------
    wav_grid : Tuple[float, float, float]
        Wavenumber grid by minimum, maximum and resolution.
    tem_grid : Tuple[int, float, float]
        Temperature grid by number of points, minimum and maximum.
    absorbers : List
        A list of absorbers.
    hitran_file : str
        Path to HITRAN file.

    Returns
    -------
    driver : Dict[str, str]
        A dictionary containing the driver file content.
    """
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

def write_rfm_atm(atm: Dict[str, np.ndarray], rundir: str=".") -> None:
    """
    Write RFM atmosphere to file.

    Parameters
    ----------
    atm : Dict[str, np.ndarray]
        A dictionary containing the atmosphere
    rundir : str
        Directory to write the file. Default is current directory.

    Returns
    -------
    None
    """
    print(f"# Creating {rundir}/rfm.atm ...")
    num_layers = atm["HGT"].shape[0]
    if not os.path.exists(f"{rundir}"):
        os.makedirs(f"{rundir}")

    with open(f"{rundir}/rfm.atm", "w") as file:
        file.write("%d\n" % num_layers)
        file.write("*HGT [km]\n")
        for i in range(num_layers):  # m -> km
            file.write("%.8g " % (atm["HGT"][i] / 1.0e3,))
        file.write("\n*PRE [mb]\n")
        for i in range(num_layers):  # pa -> mb
            file.write("%.8g " % (atm["PRE"][i] / 100.0,))
        file.write("\n*TEM [K]\n")
        for i in range(num_layers):
            file.write("%.8g " % atm["TEM"][i])
        for name, val in atm.items():
            if name in ["IDX", "HGT", "PRE", "TEM"]:
                continue
            file.write("\n*" + name + " [ppmv]\n")
            for j in range(num_layers):  # mol/mol -> ppmv
                file.write("%.8g " % (val[j] * 1.0e6,))
        file.write("\n*END")
    print(f"# {rundir}/rfm.atm written.")

def write_rfm_drv(driver: Dict[str, str], rundir: str=".") -> None:
    """
    Write RFM driver to file.

    Parameters
    ----------
    driver : Dict[str, str]
        A dictionary containing the driver file content.
    rundir : str
        Directory to write the file. Default is current directory

    Returns
    -------
    None
    """
    print(f"# Creating {rundir}/rfm.drv ...")
    if not os.path.exists(f"{rundir}"):
        os.makedirs(f"{rundir}")

    with open(f"{rundir}/rfm.drv", "w") as file:
        for sec in driver:
            if driver[sec] != None:
                file.write(sec + "\n")
                file.write(" " * 4 + driver[sec] + "\n")
    print(f"# {rundir}/rfm.drv written.")

def run_rfm(rundir: str=".") -> None:
    """
    Call to run RFM.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    pwd = os.getcwd()
    if not os.path.exists(rundir):
        os.makedirs(rundir)

    if rundir != ".":
        os.chdir(os.path.join(pwd, rundir))

    with open(f"rfm.runlog", "w") as file:
        process = subprocess.Popen(
            [f"{pwd}/rfm.release"], stdout=file, stderr=subprocess.STDOUT
        )
        process.communicate()

    #for line in iter(process.stdout.readline, b""):
        # decode the byte string and end='' to avoid double newlines
    #    print(line.decode(), end="")


def create_netcdf_input(
    fname: str,
    absorbers: List[str],
    atm: Dict[str, np.ndarray],
    wmin: float,
    wmax: float,
    wres: float,
    tnum: int,
    tmin: float,
    tmax: float,
) -> str:
    """
    Create an input file for writing kcoeff table to netCDF format

    Parameters
    ----------
    fname : str
        Name of the file.
    absorbers : list
        A list of absorbers.
    atm : Dict[str, np.ndarray]
        A dictionary containing the atmosphere.
    wmin : float
        Minimum wavenumber.
    wmax : float
        Maximum wavenumber.
    wres : float
        Wavenumber resolution.
    tnum : int
        Number of temperature points.
    tmin : float
        Minimum temperature.
    tmax : float
        Maximum temperature.

    Returns
    -------
    fname : str
        Name of the input file for netCDf
    """
    print(f"# Creating {fname}.in ...")

    with open(f"{fname}.in", "w") as file:
        file.write("# Molecular absorber\n")
        file.write("%d\n" % len(absorbers))
        file.write(" ".join(absorbers) + "\n")
        file.write("# Molecule data files\n")
        for ab in absorbers:
            file.write("%-40s\n" % (f"tab_" + ab.lower() + ".txt",))
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
    print(f"# {fname}.in written.")
    return f"{fname}.in"

def write_ktable(
    fname: str,
    absorbers: List[str],
    atm: Dict[str, np.ndarray],
    wav_grid: Tuple[float, float, float],
    tem_grid: Tuple[int, float, float],
    basedir: str=".",
) -> None:
    """
    Write kcoeff table to netCDF file.

    Parameters
    ----------
    fname : str
        Name of the file.
    absorbers : List
        A list of absorbers.
    atm : Dict[str, np.ndarray]
        A dictionary containing the atmosphere.
    wav_grid : Tuple[float, float, float]
        Wavenumber grid by minimum, maximum and resolution.
    tem_grid : Tuple[int, float, float]
        Temperature grid by number of points, minimum and maximum.

    Returns
    -------
    None
    """
    inpfile = create_netcdf_input(fname, absorbers, atm, *wav_grid, *tem_grid)

    process = subprocess.Popen(
        [f"{basedir}/kcoeff.release", "-i", inpfile, "-o", f"{fname}.nc"],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    for line in iter(process.stdout.readline, b""):
        # decode the byte string and end='' to avoid double newlines
        print(line.decode(), end="")

    process.communicate()

    pwd = os.getcwd()
    shutil.move(f"{fname}.nc", f"{basedir}/{fname}.nc")
    print(f"# {fname}.nc written.")
