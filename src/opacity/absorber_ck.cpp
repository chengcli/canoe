// C/C++
#include <algorithm>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// canoe
#include <air_parcel.hpp>
#include <constants.hpp>

// netcdf
#ifdef NETCDFOUTPUT
extern "C" {
#include <netcdf.h>
}
#endif

// climath
#include <climath/interpolation.h>

// opacity
#include "absorber_ck.hpp"

void AbsorberCK::LoadCoefficient(std::string fname, size_t bid) {
#ifdef NETCDFOUTPUT
  int fileid, dimid, varid, err;
  len_[0] = 22;  // number of pressure
  len_[1] = 27;  // number of temperature
  // len_[2] = 26; // number of bands
  len_[2] = 8;  // number of guass points

  nc_open(fname.c_str(), NC_NETCDF4, &fileid);

  axis_.resize(len_[0] + len_[1] + len_[2]);

  nc_inq_varid(fileid, "p", &varid);
  nc_get_var_double(fileid, varid, axis_.data());
  nc_inq_varid(fileid, "t", &varid);
  nc_get_var_double(fileid, varid, axis_.data() + len_[0]);
  nc_inq_varid(fileid, "samples", &varid);
  nc_get_var_double(fileid, varid, axis_.data() + len_[0] + len_[1]);

  kcoeff_.resize(len_[0] * len_[1] * len_[2]);
  nc_inq_varid(fileid, GetName().c_str(), &varid);
  size_t start[4] = {0, 0, bid, 0};
  size_t count[4] = {len_[0], len_[1], 1, len_[2]};
  nc_get_vara_double(fileid, varid, start, count, kcoeff_.data());
  nc_close(fileid);
#endif
}

Real AbsorberCK::GetAttenuation(Real g1, Real g2,
                                AirParcel const& var) const {
  // first axis is wavenumber, second is pressure, third is temperature anomaly
  Real val, coord[3] = {log(var.q[IPR]), var.q[IDN], g1};
  interpn(&val, coord, kcoeff_.data(), axis_.data(), len_, 3, 1);
  Real dens = var.q[IPR] / (Constants::kBoltz * var.q[IDN]);
  return exp(val) * dens;  // ln(m*2/kmol) -> 1/m
}
