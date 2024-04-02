#! /usr/bin/env python3
import sys, os
import numpy as np
sys.path.append("../python")
sys.path.append(".")

from canoe import def_species, load_configure
from canoe.snap import def_thermo
from canoe.athena import Mesh, ParameterInput, Outputs,MeshBlock
from typing import Tuple
from canoe.harp import radiation_band, radiation

def set_atmos_run_RT(mb: MeshBlock,
                     qNH3: float, # ppmv 
                     T0: float  # Kelvin
                     ):
    
    mb.construct_atmosphere(pin,qNH3,T0)    
    rad=mb.get_rad()
    rad.cal_radiance(mb, mb.k_st, mb.j_st)

    nb=rad.get_num_bands()
    tb=np.array([0.]*4*nb)

    for ib in range(nb):
        # print(rad.get_band(ib))
        toa=rad.get_band(ib).get_toa()[0]
        tb[ib*4:ib*4+4]=toa
    return tb


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
    rad=mb.get_rad()
    rad.cal_radiance(mb, mb.k_st, mb.j_st)

    nb=rad.get_num_bands()
    tb=np.array([0.]*4*nb)

    for ib in range(nb):
        print()
        band1=rad.get_band(ib)
        print
        toa=rad.get_band(ib).get_toa()[0]
        tb[ib*4:ib*4+4]=toa
    print(tb)

    # print(len(mesh.meshblocks()))
    # for mb in mesh.meshblocks():
    #    print(mb.block_size.nx1)

    tb=set_atmos_run_RT(mb,320,169)
    print(tb)

    exit()
    adlnTdlnP = 0.2  ##k
    pmin = 100.0e5
    pmax = 800.0e5
    mb.modify_dlnTdlnP(adlnTdlnP, pmin, pmax)

    adlnNH3dlnP = 0.25  ##ppmv  -100
    pmin = 100.0e5
    pmax = 800.0e5
    # mb.modify_dlnNH3dlnP(adlnNH3dlnP, pmin, pmax)
    xNH3=200
    T0=190
    # mb.construct_atmosphere(pin,xNH3,T0)


    # for k in range(mb.k_st, mb.k_ed + 1):
    #    for j in range(mb.j_st, mb.j_ed + 1):
    # print(str(k) + " " + str(j))
    # mb.get_rad().cal_radiance(mb, mb.k_st, mb.j_st)

    # pin.set_string("job", "problem_id", "juno_mwr_modify_Temp")
    # pin.set_string("job", "problem_id", "juno_mwr_modify_NH3")
    # print(pin.get_string("job", "problem_id"))

    # out = Outputs(mesh, pin)
    # out.make_outputs(mesh, pin)

    # run rt again
