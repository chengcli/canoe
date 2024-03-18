// C/C++
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
#include "oxygen_cia.hpp"

// netcdf
#ifdef NETCDFOUTPUT
extern "C" {
#include <netcdf.h>
}
#endif

static int const n_fd_sets = 15;
static int const n_a1dg_x3sg00_sets = 1;
static int const n_a1dg_x3sg10_sets = 1;
static int const n_ab_sets = 1;
static int const n_other_sets = 1;

void O2O2CIA::LoadCoefficient(std::string fname, int) {
  std::string line;
  std::ifstream file(fname);
  double value;

  // Fundamental band
  std::getline(file, line);
  auto sline = Vectorize<std::string>(line.c_str());
  int nwave = std::stoi(sline[3]), ntemp = n_fd_sets;

  fd_axis_.resize(ntemp + nwave);
  fd_.resize(ntemp * nwave);

  for (int i = 0; i < ntemp; ++i) {
    fd_axis_[i] = std::stod(sline[4]);
    for (int j = 0; j < nwave; ++j) {
      file >> value;
      fd_axis_[ntemp + j] = value;
      file >> value;
      fd_[i * nwave + j] = std::max(value, 0.);
    }
    std::getline(file, line);
    std::getline(file, line);
    sline = Vectorize<std::string>(line.c_str());
  }
  fd_len_[0] = ntemp;
  fd_len_[1] = nwave;

  // a1dg_x3sg00
  nwave = std::stoi(sline[3]), ntemp = n_a1dg_x3sg00_sets;

  a1dg_x3sg00_axis_.resize(ntemp + nwave);
  a1dg_x3sg00_.resize(ntemp * nwave);

  for (int i = 0; i < ntemp; ++i) {
    a1dg_x3sg00_axis_[i] = std::stod(sline[4]);
    for (int j = 0; j < nwave; ++j) {
      file >> value;
      a1dg_x3sg00_axis_[ntemp + j] = value;
      file >> value;
      a1dg_x3sg00_[i * nwave + j] = std::max(value, 0.);
    }
    std::getline(file, line);
    std::getline(file, line);
    sline = Vectorize<std::string>(line.c_str());
  }
  a1dg_x3sg00_len_[0] = ntemp;
  a1dg_x3sg00_len_[1] = nwave;

  // a1dg_x3sg10
  nwave = std::stoi(sline[3]), ntemp = n_a1dg_x3sg10_sets;

  a1dg_x3sg10_axis_.resize(ntemp + nwave);
  a1dg_x3sg10_.resize(ntemp * nwave);

  for (int i = 0; i < ntemp; ++i) {
    a1dg_x3sg10_axis_[i] = std::stod(sline[4]);
    for (int j = 0; j < nwave; ++j) {
      file >> value;
      a1dg_x3sg10_axis_[ntemp + j] = value;
      file >> value;
      a1dg_x3sg10_[i * nwave + j] = std::max(value, 0.);
    }
    std::getline(file, line);
    std::getline(file, line);
    sline = Vectorize<std::string>(line.c_str());
  }
  a1dg_x3sg10_len_[0] = ntemp;
  a1dg_x3sg10_len_[1] = nwave;

  // A band
  nwave = std::stoi(sline[3]), ntemp = n_ab_sets;

  ab_axis_.resize(ntemp + nwave);
  ab_.resize(ntemp * nwave);

  for (int i = 0; i < ntemp; ++i) {
    ab_axis_[i] = std::stod(sline[4]);
    for (int j = 0; j < nwave; ++j) {
      file >> value;
      ab_axis_[ntemp + j] = value;
      file >> value;
      ab_[i * nwave + j] = std::max(value, 0.);
      file >> value;
    }
    std::getline(file, line);
    std::getline(file, line);
    sline = Vectorize<std::string>(line.c_str());
  }
  ab_len_[0] = ntemp;
  ab_len_[1] = nwave;

  // other bands
  nwave = std::stoi(sline[3]), ntemp = n_other_sets;

  other_axis_.resize(ntemp + nwave);
  other_.resize(ntemp * nwave);

  for (int i = 0; i < ntemp; ++i) {
    other_axis_[i] = std::stod(sline[4]);
    for (int j = 0; j < nwave; ++j) {
      file >> value;
      other_axis_[ntemp + j] = value;
      file >> value;
      other_[i * nwave + j] = std::max(value, 0.);
      file >> value;
    }
    std::getline(file, line);
    std::getline(file, line);
    sline = Vectorize<std::string>(line.c_str());
  }
  other_len_[0] = ntemp;
  other_len_[1] = nwave;
}

Real O2O2CIA::getAttenuation1(Real wave, AirParcel const& qfrac) const {
  Real kk, temp = qfrac.w[IDN];

  Real coor[2] = {temp, wave};

  int iO2 = mySpeciesId(0);

  if (wave >= 1150. && wave <= 1950.) {
    interpnf(&kk, coor, fd_.data(), fd_axis_.data(), fd_len_, 2);
  } else if (wave >= 7450. && wave <= 8487.)
    interpnf(&kk, coor, a1dg_x3sg00_.data(), a1dg_x3sg00_axis_.data(),
             a1dg_x3sg00_len_, 2);
  else if (wave >= 9001. && wave <= 9997.)
    interpnf(&kk, coor, a1dg_x3sg10_.data(), a1dg_x3sg10_axis_.data(),
             a1dg_x3sg10_len_, 2);
  else if (wave >= 12600. && wave <= 13839.)
    interpnf(&kk, coor, ab_.data(), ab_axis_.data(), ab_len_, 2);
  else if (wave >= 14996. && wave <= 29790.)
    interpnf(&kk, coor, other_.data(), other_axis_.data(), other_len_, 2);
  else
    kk = 0.;

  Real num_o2 = qfrac.w[IPR] / (Constants::kBoltz * temp) * qfrac.w[iO2];

  // cm^5 molecule^{-2} * (molecule m^{-3})^2
  return 1.E-10 * num_o2 * num_o2 * kk;
}
