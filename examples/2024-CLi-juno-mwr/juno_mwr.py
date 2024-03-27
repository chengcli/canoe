#! /usr/bin/env python3
import sys, os

sys.path.append("../python")
sys.path.append(".")

from canoe import def_species, def_thermo, load_configure
from canoe.athena import Mesh, ParameterInput
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
    config = load_configure("mwr_channels.yaml")

    def_species(vapors=["H2O", "NH3"])
    # def_thermo(config["thermo"])

    pin = ParameterInput()
    pin.load_from_file("juno_mwr.inp")
    print(pin.get_real("problem", "qH2O.ppmv"))

    mesh = Mesh(pin)
    mesh.initialize(pin)

    nlyr = 100

    # pressure scale height
    pmax = 100.0e5
    pmin = 0.1e5

    P0 = 1.0e5
    T0 = 169.0

    comp = {}
    comp["NH3"] = 300.0e-6
    comp["H2O"] = 3000.0e-6

    # mb = construct_atmosphere(100, comp, T0, P0, pim=[pmin, pmax])

    # mb.add_radiation(config)
    # rad = mb.cal_radiance()

    # mb = modify_atmosphere(mb, T0, P0, plim=[pmin, pmax])
    # rad = mb.cal_radiance()
