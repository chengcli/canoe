#! /usr/bin/env python3
from canoe import *
from pyharp import radiation_band, subscribe_species
from utilities import load_file
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
    atm["CO2"] = 400.0 * ones(nlyr)
    # ppmv
    atm["H2O"] = 4000.0 * exp(-atm["HGT"] / (Hscale_km / 5.0))

    return atm


def planck(wavenumber, T):
    """
    Planck's Law as a function of wavenumber.

    Args:
    - wavenumber (float): Wavenumber (cm^-1).
    - T (float): Temperature in Kelvin.

    Returns:
    - float: Spectral radiance.
    """
    c1 = 1.191042e-5  # mW m^-2 sr^-1 cm^-4
    c2 = 1.4387752  # cm K
    return c1 * wavenumber**3 / (exp(c2 * wavenumber / T) - 1.0) * 1.0e-3


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

    atm = create_atmosphere(num_layers)
    band.set_spectral_properties(atm)

    toa = band.cal_radiance().get_toa()
    tau = band.get_tau()

    wavenumber = linspace(600, 700, int((700 - 600) / 0.01) + 1)

    print(wavenumber)
    print(len(wavenumber))
    print(toa.shape)

    print(planck(wavenumber, 300))

    # print(toa, tau.shape)
    figure(1)
    ax = axes()
    ax.plot(wavenumber, toa)
    ax.plot(wavenumber, planck(wavenumber, 300), "--")
    ax.plot(wavenumber, planck(wavenumber, 200), "--")
    ax.plot(wavenumber, planck(wavenumber, 100), "--")
    ax.set_xscale("log")
    ax.set_xlabel("Wavenumber (cm$^{-1}$)")
    ax.set_ylabel("Radiance (W m$^{-2}$ sr$^{-1}$ cm$^{-1}$)")
    show()
