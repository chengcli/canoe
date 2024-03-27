#! /usr/bin/env python3
import sys, os

sys.path.append("../python")
sys.path.append(".")

from canoe import def_species, load_configure
from canoe.snap import def_thermo
from canoe.athena import Mesh, ParameterInput, Outputs
from typing import Tuple

# from canoe.snap import air_parcel



def modify_atmosphere(
    mesh : Mesh,
    dlndlnP: float=0.,
    pmin : float=0.,
    pmax : float=0.,
    var  : str="Temp"
    ) -> Mesh:

    if dlndlnP !=0.:
        mesh.modify_Temp_profile()
    else: 
        pass
    return mesh


if __name__ == "__main__":
    pin = ParameterInput()
    pin.load_from_file("juno_mwr.inp")

    vapors = pin.get_string("species", "vapor").split(", ")
    clouds = pin.get_string("species", "cloud").split(", ")
    tracers = pin.get_string("species", "tracer").split(", ")

    def_species(vapors=vapors, clouds=clouds, tracers=tracers)
    def_thermo(pin)

    config = load_configure("juno_mwr.yaml")

    print(pin.get_real("problem", "qH2O.ppmv"))

    mesh = Mesh(pin)
    mesh.initialize(pin)

    out = Outputs(mesh, pin)
    out.make_outputs(mesh, pin)


    print(len(mesh.meshblocks()))

    mb = mesh.meshblock(0)
    # for mb in mesh.meshblocks():
    #    print(mb.block_size.nx1)

    # modify atmosphere
    # dlnTdlnP=0
    # pmin=0
    # pmax=0
    # adlnNH3dlnP=0
    # mesh=modify_atmosphere(mesh,dlnTdlnP,pmin,pmax,"Temp")
    # mesh=modify_atmosphere(mesh,adlnNH3dlnP,pmin,pmax,"NH3")



    # run rt again

