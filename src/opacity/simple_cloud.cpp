// C/C++ headers
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <cassert>  // assert

// Athena++ headers
//#include "../math_funcs.hpp"
#include "absorber.hpp"
#include "water_cloud.hpp"
#include "radiation_utils.hpp"  // getPhaseHenyeyGreenstein

SimpleCloud::SimpleCloud(RadiationBand *pband, int id, ParameterInput *pin):
      Absorber(pband, "simplecloud", id)
{
  char str[80];
  sprintf(str, "%s.%s.qext", pband->myname.c_str(), myname.c_str());
  qext = pin->GetOrAddReal("radiation", str, 1.);
  sprintf(str, "%s.%s.ssa", pband->myname.c_str(), myname.c_str());
  ww = pin->GetOrAddReal("radiation", str, 1.);
  sprintf(str, "%s.%s.asymf", pband->myname.c_str(), myname.c_str());
  gg = pin->GetOrAddReal("radiation", str, 1.);
}

// For grey cloud
Real SimpleCloud::getAttenuation(Real wave1, Real wave2, CellVariables const& var) const
{
  Real csize = 1.E0*1.e-6; // one micron size particle
  Real crho = 1.E3;
	Real qext = params_.at("qext");
//  std::cout<<"var.c[imol_]  "<<imol_<< " " <<var.c[imol_]<<std::endl;
//  return  var.c[imol_]*qext/(4./3.*csize*crho);     // -> 1/m
//put clouds 0 and precipitation 1 together
  Real totpar = var.c[imol_]+var.c[imol_+1];
  return  totpar*qext/(4./3.*csize*crho);     // -> 1/m
}

Real SimpleCloud::getSingleScateringAlbedo(Real wave1, Real wave2, CellVariables const& var) const
{
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


void SimpleCloud::getPhaseMomentum(Real *pp, Real wave1, Real wave2, CellVariables const& var, int np) const
{
    getPhaseHenyeyGreenstein(pp, 0, params.at("gg"), np);  // 0 for HENYEY_GREENSTEIN
/*
  Real totpar = var.c[imol_]+var.c[imol_+1];

  if (totpar < 1.e-20) {
    getPhaseHenyeyGreenstein(pp, 0, 0.0, np); // 0 for HENYEY_GREENSTEIN
  } else {
    getPhaseHenyeyGreenstein(pp, 0, gg, np);  // 0 for HENYEY_GREENSTEIN
  }
*/
}
