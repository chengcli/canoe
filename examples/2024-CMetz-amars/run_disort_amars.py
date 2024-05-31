#! /usr/bin/env python3
import sys, os

sys.path.append("../python")
sys.path.append(".")

from canoe import def_species, load_configure
from canoe.snap import def_thermo
from canoe.harp import radiation_band
from canoe.athena import *
from numpy import linspace, ones, exp, ndarray
from netCDF4 import Dataset
from pylab import *


# create atmosphere dictionary
def create_atmosphere(nlyr: int) -> dict:
    atm = {}

    data = Dataset("amars-main.nc", "r")
    # m
    atm["HGT"] = (data["x1f"][:-1] + data["x1f"][1:]) / 2.0

    # pa
    atm["PRE"] = data["press"][0, :, 0, 0]

    # K
    atm["TEM"] = data["temp"][0, :, 0, 0]

    # mole fraction
    atm["H2O"] = data["vapor1"][0, :, 0, 0] / 18.0 * 44

    # mole fraction
    atm["H2S"] = data["vapor2"][0, :, 0, 0] / 34.0 * 40

    atm["SO2"] = data["vapor3"][0, :, 0, 0] / 64 * 40

    atm["CO2"] = data["vapor4"][0, :, 0, 0] / 44 * 40

    return atm


if __name__ == "__main__":
    pin = ParameterInput()
    pin.load_from_file("amars.inp")

    vapors = pin.get_string("species", "vapor").split(", ")
    clouds = pin.get_string("species", "cloud").split(", ")
    # tracers = pin.get_string("species", "tracer").split(", ")

    def_species(vapors=vapors, clouds=clouds)
    def_thermo(pin)

    config = load_configure("amars.yaml")

    mesh = Mesh(pin)
    mesh.initialize(pin)

    mb = mesh.meshblock(0)
    rad = mb.get_rad()

    rad.cal_flux(mb, mb.k_st, mb.j_st, mb.i_st, mb.i_ed)
    print(ndarray(rad.fluxup))

    # band = radiation_band("B1", config, load_opacity=True)

    # aH2O = band.get_absorber_by_name("H2O")
    # aH2S = band.get_absorber_by_name("H2S")

    temp = 200.0
    xH2O = 0.003
    xH2S = 0.0003
    xSO2 = 0.0
    xCO2 = 0.0
    v1 = 0.0
    v2 = 0.0
    v3 = 0.0
    pres = 0.5e5

    # atm = [temp, xH2O, xH2S, xSO2, xCO2, v1, v2, v3, pres]
    # print(aH2O.get_attenuation(50.0, 50.0, atm))
    # print(aH2S.get_attenuation(50.0, 50.0, atm))

    # num_layers = 128
    # band.resize(num_layers)

    # atm = create_atmosphere(num_layers)
    # band.set_spectral_properties(atm)

    # dtau = band.get_dtau()
    # print(dtau.shape)
    # print(dtau)

    # band.cal_flux()
