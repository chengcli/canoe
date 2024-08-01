// C/C++ headers
#include <algorithm>
#include <cassert>  // assert
#include <cstring>
#include <iostream>
#include <stdexcept>

// canoe
#include <air_parcel.hpp>

// harp
#include <harp/radiation.hpp>

// opacity
#include "simple_cloud.hpp"

// For grey cloud
Real SimpleCloud::getAttenuation1(Real wave, AirParcel const& qfrac) const {
  AirParcel var(qfrac);
  var.ToMassConcentration();

  Real csize = 1.E0 * 1.e-6;  // one micron size particle
  Real crho = 1.E3;

  Real qext = GetPar<Real>("qext");
  // put clouds 0 and precipitation 1 together
  // Real totpar = var.c[myCloudId(0)];
  Real totpar = var.c[myCloudId(0)] + var.c[myCloudId(1)];
  if (totpar < 1e-10) {
    return 0.0;
  } else {
    //    std::cout<<"totmmr:  "<<totpar <<"  cloud: "<<var.c[myCloudId(0)]<< "
    //    rain: " <<var.c[myCloudId(1)]<<std::endl;
    return totpar * qext / (4. / 3. * csize * crho);  // -> 1/m
  }
}

Real SimpleCloud::getSingleScatteringAlbedo1(Real wave,
                                             AirParcel const& qfrac) const {
  return GetPar<Real>("ssa");
}

void SimpleCloud::getPhaseMomentum1(Real* pp, Real wave, AirParcel const& qfrac,
                                    int np) const {
  // 0 for HENYEY_GREENSTEIN
  // RadiationHelper::get_phase_momentum(pp, 0, GetPar<Real>("asymf"), np);
  RadiationHelper::get_phase_momentumDHG(
      pp, 0, GetPar<Real>("asymf"), GetPar<Real>("g1"), GetPar<Real>("g2"), np);
}
