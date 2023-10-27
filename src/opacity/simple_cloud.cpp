// C/C++ headers
#include <algorithm>
#include <cassert>  // assert
#include <cstring>
#include <iostream>
#include <stdexcept>

// canoe
#include <air_parcel.hpp>

// harp
#include <harp/radiation_utils.hpp>

// opacity
#include "water_cloud.hpp"

// For grey cloud
Real SimpleCloud::getAttenuation1(Real wave, AirParcel const& var) const {
  Real csize = 1.E0 * 1.e-6;  // one micron size particle
  Real crho = 1.E3;
  Real qext = params_.at("qext");
  //  std::cout<<"var.c[imol_]  "<<imol_<< " " <<var.c[imol_]<<std::endl;
  //  return  var.c[imol_]*qext/(4./3.*csize*crho);     // -> 1/m
  // put clouds 0 and precipitation 1 together
  Real totpar = var.c[imols_[0]] + var.c[imols_[1]];
  return totpar * qext / (4. / 3. * csize * crho);  // -> 1/m
}

Real SimpleCloud::getSingleScateringAlbedo(Real wave1, Real wave2,
                                           CellVariables const& var) const {
  return params_.at("ww");
  /*  // ssalb
    Real totpar = var.c[imol_]+var.c[imol_+1];

    if (totpar < 1.e-20) {
      return 0.0;
    } else {
    std::cout<<"var.c[imol_]  "<<ww<< " " <<var.c[imol_]<<std::endl;
    std::cout<<"var.c[imol_+1]  "<<ww<< " " <<var.c[imol_+1]<<std::endl;
      return ww;
    }
  */
}

void SimpleCloud::getPhaseMomentum(Real* pp, Real wave1, Real wave2,
                                   CellVariables const& var, int np) const {
  // 0 for HENYEY_GREENSTEIN
  get_phase_momentum(pp, 0, params.at("gg"), np);
  /*
    Real totpar = var.c[imol_]+var.c[imol_+1];

    if (totpar < 1.e-20) {
      getPhaseHenyeyGreenstein(pp, 0, 0.0, np); // 0 for HENYEY_GREENSTEIN
    } else {
      getPhaseHenyeyGreenstein(pp, 0, gg, np);  // 0 for HENYEY_GREENSTEIN
    }
  */
}
