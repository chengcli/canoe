"""
This is an example of how to run the 1D flux model.
"""

from pyharp import radiation_band, subscribe_species, init_thermo
from pyathena import nghost
from utilities import load_file
from numpy import linspace, ones, exp, array


def create_atmosphere(nlyr: int) -> dict:
    """
    Create an atmosphere dictionary containing the temperature, pressure, and
    species profiles.

    :param nlyr: number of layers
    :return: atmosphere dictionary
    """

    atm = {}
    ps_mbar = 1.0e3
    Ts_K = 300.0
    grav_ms2 = 9.8
    mu_kgmol = 2.0e-3
    Rgas_JKmol = 8.314
    Rgas_Jkg = Rgas_JKmol / mu_kgmol
    Hscale_km = Rgas_Jkg * Ts_K / grav_ms2 * 1.0e-3

    # km
    atm["HGT"] = linspace(0, 100, nlyr, endpoint=False)
    atm["PRE"] = ps_mbar * exp(-atm["HGT"] / Hscale_km)
    atm["TEM"] = Ts_K * ones(nlyr)
    # ppmv
    atm["H2O"] = 4000.0 * exp(-atm["HGT"] / (Hscale_km / 5.0))
    # ppmv
    atm["H2O(c)"] = 0.4 * (1.0 - exp(-atm["HGT"] / Hscale_km))
    # ppmv
    atm["H2"] = 1.0e6 - atm["H2O"] - atm["H2O(c)"]

    return atm


if __name__ == "__main__":
    # Enroll species internally used by pyharp
    subscribe_species(
        {
            "vapor": ["H2O"],
            "cloud": ["H2O(c)", "H2O(p)"],
        }
    )

    # number of atmopheric layers
    num_layers = 40 + 2 * nghost()

    # create an atmosphere dictionary
    atm = create_atmosphere(num_layers)

    info = load_file("hywater.yaml")

    init_thermo(info["thermodynamics"])

    band_names = ["ir", "vis"]
    bands = [radiation_band(name, info) for name in band_names]

    bflxup, bflxdn = [], []

    for band in bands:
        print(band)
        band.resize(num_layers, nstr=8)
        band.set_spectral_properties(atm)
        band.cal_flux()

        bflxup.append(array(band.bflxup))
        bflxdn.append(array(band.bflxdn))

    # print(bflxup, bflxdn)
