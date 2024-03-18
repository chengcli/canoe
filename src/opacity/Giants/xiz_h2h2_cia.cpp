// C/C++
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

// application
#include <application/exceptions.hpp>

// canoe
#include <air_parcel.hpp>
#include <constants.hpp>

// climath
#include <climath/core.h>
#include <climath/interpolation.h>

// utils
#include <utils/fileio.hpp>

// opacity
#include "hydrogen_cia.hpp"

XizH2H2CIA::XizH2H2CIA() : Absorber("H2-H2-CIA") { SetPar("xHe", 0.); }

void XizH2H2CIA::LoadCoefficient(std::string fname, int) {
  if (!FileExists(fname)) {
    throw NotFoundError("XizH2H2CIA", fname);
  }

  len_[1] = GetNumCols(fname) - 1;
  len_[0] = GetNumRows(fname) - 1;

  std::ifstream infile(fname.c_str(), std::ios::in);
  axis_.resize(len_[0] + len_[1]);
  kcoeff_.resize(len_[0] * len_[1]);

  Real junk;
  if (infile.is_open()) {
    infile >> junk;
    for (int j = 0; j < len_[1]; j++) {
      infile >> axis_[len_[0] + j];
    }
    for (int k = 0; k < len_[0]; k++) {
      infile >> axis_[k];
      for (int j = 0; j < len_[1]; j++) infile >> kcoeff_[k * len_[1] + j];
    }
    infile.close();
  } else {
    throw RuntimeError("XizH2H2CIA::LoadCoefficient",
                       "Cannot open file: " + fname);
  }
}

Real XizH2H2CIA::getAttenuation1(Real wave, AirParcel const& var) const {
  // first axis is wavenumber, second is temperature
  Real val, coord[2] = {wave, var.w[IDN]};

  interpn(&val, coord, kcoeff_.data(), axis_.data(), len_, 2, 1);

  Real x0 = 1.;
  if (mySpeciesId(0) == 0) {
    for (int n = 1; n <= NVAPOR; ++n) x0 -= var.w[n];
  } else {
    x0 = var.w[mySpeciesId(0)];
  }

  Real xHe = GetPar<Real>("xHe");
  Real amagat =
      x0 * var.w[IPR] / (Constants::kBoltz * var.w[IDN] * Constants::Lo);
  Real amagat_H2 = amagat * (1. - xHe);

  return 100. * exp(-val) * amagat_H2 * amagat_H2;  // 1/cm -> 1/m
}
