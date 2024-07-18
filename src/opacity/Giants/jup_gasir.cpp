// C/C++ headers
#include <algorithm>
#include <cassert>  // assert
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>

// athena
#include <athena/athena_arrays.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <air_parcel.hpp>
#include <constants.hpp>
#include <impl.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// opacity
#include "grey_gas.hpp"

// xiz semigrey for Jupiter
Real JupGasir::GetAttenuation(Real wave1, Real wave2,
                                    AirParcel const &var) const {
  Real scale = GetPar<Real>("scale");
  Real p = var.w[IPR];
  Real T = var.w[IDN];
  auto pthermo = Thermodynamics::GetInstance();
  Real dens = p / (pthermo->GetRd() * T);  // kg/m^3

  Real  jstrat = 8.e-4*pow(p, -0.5); //IR opacity from hydrocarbons and haze
  Real  cia = 2.e-8*p;

  return scale * dens * (cia + jstrat);     // -> 1/m
}
