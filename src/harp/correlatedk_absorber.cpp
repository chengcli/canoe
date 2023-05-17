// C/C++ headers
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cstring>

#include <climath/interpolation.h>
#include <snap/mesh/cell_variables.hpp>
#include <snap/constants.hpp>
//#include <snap/thermodynamics/thermodynamics.hpp>
#include "correlatedk_absorber.hpp"

// External library headers
#if NETCDFOUTPUT
  #include <netcdf.h>
#endif

void CorrelatedKAbsorber::loadCoefficient(std::string fname, int bid)
{
#if NETCDFOUTPUT
  int fileid, dimid, varid, err;
  len_[0] = 22; //number of pressure
  len_[1] = 27;  //number of temperature
  //len_[2] = 26; // number of bands
  len_[2] = 8; // number of guass points

  nc_open(fname.c_str(), NC_NETCDF4, &fileid);

  axis_.resize(len_[0] + len_[1] + len_[2]);

  nc_inq_varid(fileid, "p", &varid);
  nc_get_var_double(fileid, varid, axis_.data());
  nc_inq_varid(fileid, "t", &varid);
  nc_get_var_double(fileid, varid, axis_.data() + len_[0]);
  nc_inq_varid(fileid, "samples", &varid);
  nc_get_var_double(fileid, varid, axis_.data() + len_[0] + len_[1]);

  kcoeff_.resize(len_[0]*len_[1]*len_[2]);
  nc_inq_varid(fileid, name_.c_str(), &varid);
  size_t start[4] = {0, 0, (size_t)bid, 0};
  size_t count[4]  = {(size_t)len_[0], (size_t)len_[1], 1, (size_t)len_[2]};
  nc_get_vara_double(fileid, varid, start, count, kcoeff_.data());
  nc_close(fileid);

//  std::cout << " len: " << len_[0] <<" "<< len_[1] <<" "<< len_[2] <<" "<< len_[3] <<std::endl;
//  std::cout << " axislen: " <<axis_[len_[0]]<<std::endl;
//  for (int i = 1; i < 20; ++i) {
//	std::cout << " read: "<<kcoeff_[i*100] <<" "<<std::endl;
//  }
#endif
}

Real CorrelatedKAbsorber::getAttenuation(Real g1, Real g2, CellVariables const& var) const
{
#if NETCDFOUTPUT == OFF
  std::stringstream msg;
  msg << "### FATAL ERROR in function CorrelatedKAbsorber::getAttenuation"
      << std::endl << "NETCDF library is needed for Correlated-K absorber";
  ATHENA_ERROR(msg);
#endif
  // first axis is wavenumber, second is pressure, third is temperature anomaly
  Real val, coord[3] = {log(var.q[IPR]), var.q[IDN], g1};
  interpn(&val, coord, kcoeff_.data(), axis_.data(), len_, 3, 1);
/*
//  std::cout << " len: " << len_[0] <<" "<< len_[1] <<" "<< len_[2] <<" "<< len_[3] <<std::endl;
  std::cout << " axis: " << axis_[len_[0] + len_[1] + mw] <<" "<< axis_[len_[0] + len_[1] + len_[2] + mg] <<" "<< mw <<" "<< mg <<std::endl;
  std::cout << " inter: " << coord[0] <<" "<< coord[1] <<" "<< coord[2] <<" "<< coord[3] <<std::endl;
  for (int i = 0; i < 10; ++i) { 
  Real yy = float(i+1)*200.;
  coord[0] = log(1.E5);
  coord[1] = yy;
  interpn(&val, coord, kcoeff_.data(), axis_.data(), len_, 4);
	std::cout << " test: " << yy <<" "<<exp(val) <<" "<<std::endl;
  }
*/
  Real dens = var.q[IPR]/(Constants::kBoltz*var.q[IDN]);
//	std::cout << " test: " << dens <<" "<<exp(val) <<" "<< coord[0] <<" "<< coord[1] <<" "<< coord[2] <<" "<< coord[3] <<std::endl;
/*
  Real x0 = 1.;
  if (imol_ == 0) 
    for (int n = 1; n <= NVAPOR; ++n) x0 -= prim[n];
  else
    x0 = prim[imol_];
  return 1.E-3*exp(val)*dens*x0*mixr_; // ln(m*2/kmol) -> 1/m
*/
  return exp(val)*dens; // ln(m*2/kmol) -> 1/m
}
