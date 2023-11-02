// C/C++
#include <algorithm>
#include <cassert>  // assert
#include <cstring>
#include <iostream>
#include <stdexcept>

// canoe
#include <air_parcel.hpp>

// climath
#include <climath/interpolation.h>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// harp
#include <harp/radiation_utils.hpp>  // get_phase_momentum

// opacity
#include "water_cloud.hpp"

// coefficient from fu and liou

// wavenumber range
Real wband[19] = {
    0,    280,  400,  540,  670,  800,  980,  1100,  1250, 1400,
    1700, 1900, 2200, 2857, 4000, 5263, 7692, 14286, 50000};  // cm-1

// ap for tau

Real ap[18][3] = {{-0.67163E-03, 0.33056E+01, 0.0},
                  {0.25307E-03, 0.32490E+01, 0.0},
                  {-0.75524E-03, 0.33083E+01, 0.0},
                  {-0.20332E-02, 0.33865E+01, 0.0},
                  {0.40939E-02, 0.29870E+01, 0.0},
                  {-0.27583E-02, 0.34436E+01, 0.0},
                  {-7.770e-3, 3.734, 11.85},
                  {-8.088e-3, 3.717, 17.17},
                  {-8.441e-3, 3.715, 19.48},
                  {-9.061e-3, 3.741, 26.48},
                  //    c-- changed by Z.F. for the windows domain in longwave
                  //    spectral.^M 1            0.160239,  0.495375,  -4.38738,
                  //    1            0.165637, -0.438836,   1.54020,
                  //    1            0.172217,  -1.49513,  10.56623,
                  {-9.609e-3, 3.768, 34.11},
                  {-1.153e-2, 4.109, 17.32},
                  {-8.294e-3, 3.925, 1.315},
                  {-1.026e-2, 4.105, 16.36},
                  {-1.151e-2, 4.182, 31.13},
                  {-1.704e-2, 4.830, 16.27},
                  {-1.741e-2, 5.541, -58.42},
                  {-7.752e-3, 4.624, -42.01}};

// bp for ssalb

Real bp[18][4] = {{-0.14661E-06, 0.79495E-07, -0.10422E-09, 0.40232E-12},
                  {-0.15417E-05, 0.11489E-04, -0.77147E-08, 0.22160E-10},
                  {-0.13287E-02, 0.91493E-03, -0.39410E-05, 0.12610E-07},
                  {-0.21311E-02, 0.22827E-02, -0.13400E-04, 0.42169E-07},
                  {0.22764E+00, 0.21902E-02, -0.16743E-04, 0.53032E-07},
                  {0.59555E-01, 0.73777E-02, -0.66056E-04, 0.21750E-06},
                  {.19960E+00, .37800E-02, -.14910E-04, .00000E+00},
                  {.30140E+00, .26390E-02, -.11160E-04, .00000E+00},
                  {.39080E+00, .12720E-02, -.55640E-05, .00000E+00},
                  {.31050E+00, .26030E-02, -.11390E-04, .00000E+00},
                  {0.236894, 2.10402E-03, -3.72955E-06, 0.0},
                  {0.315225, 9.38232E-04, 1.50649E-06, 0.0},
                  {0.605243, -3.92611E-03, 2.12776E-05, 0.0},
                  //    {.20370E+00,  .42470E-02, -.18100E-04,  .00000E+00},
                  //   {.23070E+00,  .38300E-02, -.16160E-04,  .00000E+00},
                  //   {.56310E+00, -.14340E-02,  .62980E-05,  .00000E+00},
                  {.52070E+00, -.97780E-03, .37250E-05, .00000E+00},
                  {.32540E+00, .34340E-02, -.30810E-04, .91430E-07},
                  {.10280E+00, .50190E-02, -.20240E-04, .00000E+00},
                  {.39640E+00, -.31550E-02, .64170E-04, -.29790E-06},
                  {.80790E+00, -.70040E-02, .52090E-04, -.14250E-06}};

// cpir for gg   ** not right for shortwave band (iband =1 ~ 6)

Real cpir[18][4] = {
    {.79550, 2.524e-3, -1.022e-5, 0.000e+0},
    {.79550, 2.524e-3, -1.022e-5, 0.000e+0},
    {.79550, 2.524e-3, -1.022e-5, 0.000e+0},
    {.79550, 2.524e-3, -1.022e-5, 0.000e+0},
    {.79550, 2.524e-3, -1.022e-5, 0.000e+0},
    {.79550, 2.524e-3, -1.022e-5, 0.000e+0},
    {.79550, 2.524e-3, -1.022e-5, 0.000e+0},
    {.86010, 1.599e-3, -6.465e-6, 0.000e+0},
    {.89150, 1.060e-3, -4.171e-6, 0.000e+0},
    {.87650, 1.198e-3, -4.485e-6, 0.000e+0},
    {0.884846, 7.52769E-05, 4.57733E-06, 0.0},
    {0.901327, 2.03758E-04, 2.95010E-06, 0.0},
    {0.873900, 1.45318E-03, -6.30462E-06, 0.0},
    //             .88150,     9.858e-4,    -3.116e-6,     0.000e+0,
    //              .91670,     5.499e-4,    -1.507e-6,     0.000e+0,
    //              .90920,     9.295e-4,    -3.877e-6,     0.000e+0,
    {.84540, 1.429e-3, -5.859e-6, 0.000e+0},
    {.76780, 2.571e-3, -1.041e-5, 0.000e+0},
    {.72900, 2.132e-3, -5.584e-6, 0.000e+0},
    {.70240, 4.581e-3, -3.054e-5, 6.684e-8},
    {.22920, 1.724e-2, -1.573e-4, 4.995e-7}};

// TODO(cli): number density and mass density are incorrect
Real FuWaterIceCloud::getAttenuation1(Real wave, AirParcel const& qfrac) const {
  auto pthermo = Thermodynamics::GetInstance();

  // from fu & liou code
  int iband = locate(wband, wave, 19) + 1;
  assert(iband >= 1 && iband <= 18);
  iband = 18 - iband;

  Real result = 0.;
  Real fw1, fw2, pde, T;
  T = qfrac.w[IDN] - 273.15;  // K to degree C
  if (T > -60.0) {
    pde = 326.3 + 12.42 * T + 0.197 * T * T +
          0.0012 * T * T * T;  // effective size of ice cloud ( um )
  } else {
    pde = 30.;
  }
  fw1 = pde;
  fw2 = fw1 * pde;
  Real dens =
      qfrac.w[IPR] / (pthermo->GetRd() * qfrac.w[IDN] * pthermo->RovRd(qfrac));

  result = ap[iband - 1][0] + ap[iband - 1][1] / fw1 + ap[iband - 1][2] / fw2;

  return 1e3 * dens * qfrac.c[GetCloudIndex(0)] * result;
}

Real FuWaterIceCloud::getSingleScatteringAlbedo1(Real wave,
                                                 AirParcel const& qfrac) const {
  // from fu & liou code
  int iband = locate(wband, wave, 19) + 1;
  assert(iband >= 1 && iband <= 18);
  iband = 18 - iband;

  Real wi = 0.;
  Real fw1, fw2, fw3, pde, T;
  T = qfrac.w[IDN] - 273.15;  // K to degree C
  if (T > -60.0) {
    pde = 326.3 + 12.42 * T + 0.197 * T * T + 0.0012 * T * T * T;
  } else {
    pde = 30.;
  }
  fw1 = pde;
  fw2 = fw1 * pde;
  fw3 = fw2 * pde;

  wi = 1. - (bp[iband - 1][0] + bp[iband - 1][1] * fw1 +
             bp[iband - 1][2] * fw2 + bp[iband - 1][3] * fw3);

  // if (prim[imols_[0]]<2e-19){
  //   return 0.0;
  // }else{
  return wi;
  //}
}

void FuWaterIceCloud::getPhaseMomentum1(Real* pp, Real wave,
                                        AirParcel const& qfrac, int np) const {
  // from fu & liou code
  int iband = locate(wband, wave, 19) + 1;
  assert(iband >= 1 && iband <= 18);
  iband = 18 - iband;

  Real fw1, fw2, fw3, pde, T;
  T = qfrac.w[IDN] - 273.15;  // K to degree C
  if (T > -60.0 && T < 0.) {
    pde = 326.3 + 12.42 * T + 0.197 * T * T + 0.0012 * T * T * T;
  } else {
    pde = 30.;
  }

  fw1 = pde;
  fw2 = fw1 * pde;
  fw3 = fw2 * pde;

  Real gg = 0.0;
  gg = cpir[iband - 1][0] + cpir[iband - 1][1] * fw1 +
       cpir[iband - 1][2] * fw2 + cpir[iband - 1][3] * fw3;

  // std::cout<<"pde "<<pde<<" "<<T<<" "<<gg<<"ice"<<std::endl;
  // if (prim[imols_[0]]<2e-19){
  //   GetMom(0, 0.0, np, pp); // 0 for HENYEY_GREENSTEIN
  // }else{
  RadiationHelper::get_phase_momentum(pp, 0, gg, np);
  //}
}
