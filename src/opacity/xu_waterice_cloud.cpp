// C/C++ headers
#include <algorithm>
#include <cassert>  // assert
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>

// canoe
#include <air_parcel.hpp>

// climath
#include <climath/interpolation.h>

// utils
#include <utils/vectorize.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// harp
#include <harp/radiation.hpp>  // GetPhaseMomentum

// opacity
#include "water_cloud.hpp"

void XuWaterIceCloud::LoadCoefficient(std::string fname, size_t bid) {
  std::string line;
  std::ifstream file(fname.c_str(), std::ios::in);

  int nwave = 396, nsize = 189;

  axis_.resize(nwave + nsize);
  ssalb_.resize(nwave * nsize);
  gg_.resize(nwave * nsize);

  Real value;
  for (int i = 0; i < nwave; ++i) {
    for (int j = 0; j < nsize; ++j) {
      std::getline(file, line);

      auto sline = Vectorize<std::string>(line.c_str());
      if (i == 0) axis_[nwave + j] = atof(sline[1].c_str());
      if (j == 0) axis_[i] = atof(sline[0].c_str());
      ssalb_[i * nsize + j] = atof(sline[5].c_str());
      gg_[i * nsize + j] = atof(sline[6].c_str());
      /*
      if (i ==0) axis_[nwave+j] = std::stod(sline[1]);
      if (j ==0) axis_[i]  = std::stod(sline[0]);
      ssalb_[i*nsize+j] = std::stod(sline[5]);
      gg_[i*nsize+j] = std::stod(sline[6]);
      */
    }
  }
  len_[0] = nwave;  // number of wavelength
  len_[1] = nsize;  // number of size
}

Real XuWaterIceCloud::getAttenuation1(Real wave, AirParcel const& air) const {
  Real T1, re, result;
  // calculate ice cloud size
  T1 = air.w[IDN] - 273.15;  // K to C
  T1 = std::min(-25., std::max(-50., T1));
  re = 326.3 + 12.42 * T1 + 0.197 * T1 * T1 + 0.0012 * T1 * T1 * T1;
  Real dens = air.w[IPR] * 29e-3 / (Constants::Rgas * air.w[IDN]);

  result = (0.003448 + 2.431 / re) * dens * air.c[GetCloudIndex(0)] * 1e9;

  return result;
}

Real XuWaterIceCloud::getSingleScatteringAlbedo1(Real wave,
                                                 AirParcel const& air) const {
  Real T1, re, kk;
  // calculate ice cloud size
  T1 = air.w[IDN] - 273.15;  // K to C
  T1 = std::min(-25., std::max(-50., T1));
  re = 326.3 + 12.42 * T1 + 0.197 * T1 * T1 + 0.0012 * T1 * T1 * T1;

  Real ww = 10000. / wave;  // wavelength(um) to wavenumber (cm-1)
  Real coor[2] = {ww, re};

  interpnf(&kk, coor, ssalb_.data(), axis_.data(), len_, 2);

  // if (air.w[imols_[0]]<2e-19){
  //   return 0.0;
  // }else{
  return kk;
  //}
}

void XuWaterIceCloud::getPhaseMomentum1(Real* pp, Real wave,
                                        AirParcel const& air, int np) const {
  Real T1, re, gg;
  // calculate ice cloud size
  T1 = air.w[IDN] - 273.15;  // K to C
  T1 = std::min(-25., std::max(-50., T1));
  re = 326.3 + 12.42 * T1 + 0.197 * T1 * T1 + 0.0012 * T1 * T1 * T1;

  Real ww = 10000. / wave;  // wavelength(um) to wavenumber (cm-1)
  Real coor[2] = {ww, re};

  interpnf(&gg, coor, gg_.data(), axis_.data(), len_, 2);

  // if (air.w[imols_[0]]<2e-19){
  //   get_phase_momentum(pp, 0.0, np); // 0 for HENYEY_GREENSTEIN
  // }else{
  RadiationHelper::get_phase_momentum(pp, 0, gg, np);
  //}
}
