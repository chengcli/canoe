#! /usr/bin/env python3
import sys, os

sys.path.append("../python")
sys.path.append(".")

from canoe import def_species, load_configure
from canoe.snap import def_thermo
from canoe.athena import Mesh, ParameterInput, Outputs
from typing import Tuple

# from canoe.snap import air_parcel


def construct_atmosphere(
    nlyr: int = 100,
    comp: dict = {},
    plim: Tuple[float, float] = [0.1e5, 100.0e5],
    zlim: Tuple[float, float] = [0.0, 10.0e3],
    T0: float = 169.0,
    P0: float = 1.0e5,
    Tmin: float = 100.0,
) -> Mesh:
    mb = Mesh(nlyr)

    air = air_parcel(comp)
    dlnp = (log(plim[1]) - log(plim[0])) / nlyr

    mb.set_layer(0, air)
    for i in range(1, nlyr):
        atm_extrapolate(air, -dlnp / 2.0, method="dry_adiabat")
        mb.set_layer(i, air)

    return mb


def modify_atmosphere(
    mb: Mesh,
    comp: dict = {},
    plim: Tuple[float, float] = [0.1e5, 100.0e5],
    zlim: Tuple[float, float] = [0.0, 10.0e3],
    T0: float = 169.0,
    P0: float = 1.0e5,
) -> Mesh:
    pass


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
