#! /usr/bin/env python3
import sys, os

sys.path.append("../python")
sys.path.append(".")

import canoe
from canoe import load_configure
from canoe.harp import radiation_band
from atm_profile_utils import read_atm_profile
import numpy as np

def plot_opacity():
    config = load_configure("amarsw-rt.yaml")
    band = radiation_band("B1", config, load_opacity=True)
    aH2O = band.get_absorber_by_name("H2O")
    
    atm = read_atm_profile("amarsw.atm")
    ilyr = 0
    atm = [atm['TEM'][ilyr], atm['H2O'][ilyr], 0., 0., 0., atm['PRE'][ilyr]]
    print(atm)
    print(aH2O.get_attenuation(50.0, 50.0, atm))

if __name__ == "__main__":
    mesh, inp = canoe.start_with_input("amarsw.inp")

    mb = mesh.meshblock(0)
    rad = mb.get_rad()

    atm = read_atm_profile("amarsw.atm")
    mb.modify_atm(atm)

    rad.cal_flux(mb)
    x = np.array(rad.fluxup)
    y = np.array(rad.fluxdn)
    print(x.sum(axis = (0,)))
    print(y.sum(axis = (0,)))
