#! /usr/bin/env python3
import sys, os

sys.path.append("../python")
sys.path.append(".")

from canoe import def_species, load_configure
from canoe.snap import def_thermo
from canoe.athena import Mesh, ParameterInput, Outputs
from typing import Tuple
from canoe.harp import radiation_band, radiation

# from canoe.snap import air_parcel

if __name__ == "__main__":
    pin = ParameterInput()
    pin.load_from_file("juno_mwr.inp")

    vapors = pin.get_string("species", "vapor").split(", ")
    clouds = pin.get_string("species", "cloud").split(", ")
    tracers = pin.get_string("species", "tracer").split(", ")

    def_species(vapors=vapors, clouds=clouds, tracers=tracers)
    def_thermo(pin)

    config = load_configure("juno_mwr.yaml")

    # print(pin.get_real("problem", "qH2O.ppmv"))

    mesh = Mesh(pin)
    mesh.initialize(pin)

    # out = Outputs(mesh, pin)
    # out.make_outputs(mesh, pin)

    # exit()
    mb = mesh.meshblock(0)
    # print(len(mesh.meshblocks()))
    # for mb in mesh.meshblocks():
    #    print(mb.block_size.nx1)

    adlnTdlnP = 0.2  ##k
    pmin = 100.0e5
    pmax = 800.0e5
    mb.modify_dlnTdlnP(adlnTdlnP, pmin, pmax)

    adlnNH3dlnP = 0.25  ##ppmv  -100
    pmin = 100.0e5
    pmax = 800.0e5
    # mb.modify_dlnNH3dlnP(adlnNH3dlnP, pmin, pmax)

    prad = radiation(mb, pin)

    for k in range(mb.k_st, mb.k_ed + 1):
        for j in range(mb.j_st, mb.j_ed + 1):
            print(str(k) + " " + str(j))
            prad.cal_radiance(mb, k, j)

    pin.set_string("job", "problem_id", "juno_mwr_modify_Temp")
    # pin.set_string("job", "problem_id", "juno_mwr_modify_NH3")
    print(pin.get_string("job", "problem_id"))

    out = Outputs(mesh, pin)
    out.make_outputs(mesh, pin)

    # run rt again
