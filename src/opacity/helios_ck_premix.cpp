// C/C++
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

// canoe
#include <air_parcel.hpp>
#include <constants.hpp>

// climath
#include <climath/interpolation.h>

// opacity
#include "absorber_ck.hpp"

void HeliosCKPremix::LoadCoefficient(std::string fname, size_t bid) {
  std::ifstream file(fname);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file: " + fname);
  }

  // skip the first line
  std::string junk_line;
  std::getline(file, junk_line);

  size_t num_bands;

  // temperature, pressure, band, g-points
  file >> len_[0] >> len_[1] >> num_bands >> len_[2];

  if (bid >= num_bands) {
    throw std::runtime_error("Band index out of range: " + std::to_string(bid));
  }

  axis_.resize(len_[0] + len_[1] + len_[2]);
  kcoeff_.resize(len_[0] * len_[1] * len_[2]);

  // temperature grid
  for (int i = 0; i < len_[0]; ++i) {
    file >> axis_[i];
  }

  // pressure grid
  for (int j = 0; j < len_[1]; ++j) {
    Real pres;
    file >> pres;
    axis_[len_[0] + j] = log(pres);
  }

  Real wmin, wmax;
  file >> wmin;
  // band wavelengths
  for (int b = 0; b < num_bands; ++b) {
    if (b == bid) {
      file >> wmax;
      break;
    } else {
      file >> wmin;
    }
  }

  // skip unimportant wavelengths
  Real dummy;
  for (int i = 0; i < num_bands - bid; ++i) {
    file >> dummy;
  }

  // g-points and weights
  for (int g = 0; g < len_[2]; ++g) {
    Real gpoint;
    file >> gpoint >> dummy;
    axis_[len_[0] + len_[1] + g] = wmin + (wmax - wmin) * gpoint;
  }

  int n = 0;
  for (int i = 0; i < len_[0]; ++i)
    for (int j = 0; j < len_[1]; ++j)
      for (int b = 0; b < num_bands; ++b) {
        if (b == bid) {
          for (int g = 0; g < len_[2]; ++g, ++n) {
            file >> kcoeff_[n];
            kcoeff_[n] = log(std::max(kcoeff_[n], 1.0e-99));
          }
        } else {
          for (int g = 0; g < len_[2]; ++g, ++n) {
            file >> dummy;
          }
        }
      }

  file.close();
}

Real HeliosCKPremix::GetAttenuation(Real g1, Real g2,
                                    AirParcel const& var) const {
  // temperature, log-pressure, wave-scaled g-point
  Real val, coord[3] = {var.q[IDN], log(var.q[IPR]), g1};
  interpn(&val, coord, kcoeff_.data(), axis_.data(), len_, 3, 1);
  Real dens = var.q[IPR] / (Constants::kBoltz * var.q[IDN]);
  return exp(val) * dens;  // ln(cm^2 / molecule) -> 1/m
}
