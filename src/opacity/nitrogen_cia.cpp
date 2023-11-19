// C/C++
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

// canoe
#include <air_parcel.hpp>
#include <configure.hpp>
#include <constants.hpp>

// climath
#include <climath/interpolation.h>

// utils
#include <utils/vectorize.hpp>

// opacity
#include "nitrogen_cia.hpp"

// netcdf
#ifdef NETCDFOUTPUT
extern "C" {
#include <netcdf.h>
}
#endif

static int const n_rt_sets = 10;
static int const n_fd1_sets = 5;
static int const n_fd2_sets = 5;

void N2N2CIA::LoadCoefficient(std::string fname, size_t bid) {
  std::string line;
  std::ifstream file(fname);
  double value;

  // Roto-translational band
  std::getline(file, line);
  auto sline = Vectorize<std::string>(line.c_str());
  int nwave = std::stoi(sline[3]), ntemp = n_rt_sets;

  rt_axis_.resize(ntemp + nwave);
  rt_.resize(ntemp * nwave);

  for (int i = 0; i < ntemp; ++i) {
    rt_axis_[i] = std::stod(sline[4]);
    for (int j = 0; j < nwave; ++j) {
      file >> value;
      rt_axis_[ntemp + j] = value;
      file >> value;
      rt_[i * nwave + j] = std::max(value, 0.);
    }
    std::getline(file, line);
    std::getline(file, line);
    sline = Vectorize<std::string>(line.c_str());
  }
  rt_len_[0] = ntemp;
  rt_len_[1] = nwave;

  // Fundamental band
  nwave = std::stoi(sline[3]), ntemp = n_fd1_sets;

  fd1_axis_.resize(ntemp + nwave);
  fd1_.resize(ntemp * nwave);

  for (int i = 0; i < ntemp; ++i) {
    fd1_axis_[i] = std::stod(sline[4]);
    for (int j = 0; j < nwave; ++j) {
      file >> value;
      fd1_axis_[ntemp + j] = value;
      file >> value;
      fd1_[i * nwave + j] = std::max(value, 0.);
    }
    std::getline(file, line);
    std::getline(file, line);
    sline = Vectorize<std::string>(line.c_str());
  }
  fd1_len_[0] = ntemp;
  fd1_len_[1] = nwave;

  // Fundamental band
  nwave = std::stoi(sline[3]), ntemp = n_fd2_sets;

  fd2_axis_.resize(ntemp + nwave);
  fd2_.resize(ntemp * nwave);

  for (int i = 0; i < ntemp; ++i) {
    fd2_axis_[i] = std::stod(sline[4]);
    for (int j = 0; j < nwave; ++j) {
      if (i == 2) {
        file >> value;
        fd2_axis_[ntemp + j] = value;
        file >> value;
        fd2_[i * nwave + j] = std::max(value, 0.);
      } else {
        file >> value;
        fd2_axis_[ntemp + j] = value;
        file >> value;
        fd2_[i * nwave + j] = std::max(value, 0.);
        file >> value;
      }
    }
    std::getline(file, line);
    std::getline(file, line);
    sline = Vectorize<std::string>(line.c_str());
  }
  fd2_len_[0] = ntemp;
  fd2_len_[1] = nwave;
}

Real N2N2CIA::getAttenuation1(Real wave, AirParcel const& qfrac) const {
  Real kk, temp = qfrac.w[IDN];

  Real coor[2] = {temp, wave};

  int iN2 = mySpeciesId(0);

  if (wave >= 0.02 && wave <= 554.) {
    interpnf(&kk, coor, rt_.data(), rt_axis_.data(), rt_len_, 2);
  } else if (wave >= 2000. && wave <= 2698.) {
    interpnf(&kk, coor, fd2_.data(), fd2_axis_.data(), fd2_len_, 2);
  } else if (wave >= 1850. && wave <= 3000.)
    interpnf(&kk, coor, fd1_.data(), fd1_axis_.data(), fd1_len_, 2);
  else
    kk = 0.;

  Real num_n2 = qfrac.w[IPR] / (Constants::kBoltz * temp) * qfrac.w[iN2];

  // cm^5 molecule^{-2} * (molecule m^{-3})^2
  return 1.E-10 * num_n2 * num_n2 * kk;
}
