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
Real JupGasv::GetAttenuation(Real wave1, Real wave2,
                                    AirParcel const &var) const {
  Real scale = GetPar<Real>("scale");
  Real p = var.w[IPR];
  Real T = var.w[IDN];
  auto pthermo = Thermodynamics::GetInstance();
  Real dens = p / (pthermo->GetRd() * T);  // kg/m^3

//this one is a good haze
  //Real result = 1.e-6*pow(p,0.5)+1.e-3*pow(p/1.e3, -2.); //visible opacity with Jupiter haze
  //Real strongch4 = 1.e-2*pow(p, -0.5); //visible opacity with Jupiter haze
  Real strongch4 = 5.e-3*pow(p, -0.5); //visible opacity with Jupiter haze
  Real weakch4 = 0.;//1.e-3; //visible opacity with Jupiter haze
  //Real weakch4 = 1.e-3; //visible opacity with Jupiter haze

  //std::cout<<"scale=  " <<scale<<"  pres=  "<<p<< "  dens= " <<dens<<std::endl;

  return scale * dens * (strongch4+weakch4);     // -> 1/m
}
