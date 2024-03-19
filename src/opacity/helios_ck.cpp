// C/C++
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

// canoe
#include <air_parcel.hpp>
#include <constants.hpp>

// application
#include <application/exceptions.hpp>

// climath
#include <climath/interpolation.h>

// harp
#include <harp/spectral_grid.hpp>

// opacity
#include "absorber_ck.hpp"

HeliosCK::HeliosCK(std::string name) : AbsorberCK(name) {}

void HeliosCK::ModifySpectralGrid(std::vector<SpectralBin>& spec) const {
  std::cout << "I'm modifying the spectral grid" << std::endl;
  spec.resize(weights_.size());

  for (size_t i = 0; i < weights_.size(); ++i) {
    spec[i].wav1 = axis_[len_[0] + len_[1] + i];
    spec[i].wav2 = axis_[len_[0] + len_[1] + i];
    spec[i].wght = weights_[i];
    std::cout << spec[i].wav1 << " " << spec[i].wav2 << " " << spec[i].wght
              << std::endl;
  }
}

void HeliosCK::LoadCoefficient(std::string fname, int bid) {
  std::ifstream file(fname);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file: " + fname);
  }

  // over ride band id if it is set
  if (HasPar("band_id")) {
    bid = static_cast<int>(GetPar<Real>("band_id"));
  }

  // skip the first line
  std::string junk_line;
  std::getline(file, junk_line);

  int num_bands;

  // temperature, pressure, band, g-points
  file >> len_[0] >> len_[1] >> num_bands >> len_[2];

  if (bid >= num_bands || bid < 0) {
    throw RuntimeError("HeliosCK::LoadCoefficient",
                       "Band index out of range: " + std::to_string(bid));
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
  for (int i = 1; i < num_bands - bid; ++i) {
    file >> dummy;
  }

  // g-points and weights
  weights_.resize(len_[2]);
  for (int g = 0; g < len_[2]; ++g) {
    Real gpoint;
    file >> gpoint;
    axis_[len_[0] + len_[1] + g] = wmin + (wmax - wmin) * gpoint;
  }

  for (int g = 0; g < len_[2]; ++g) {
    file >> weights_[g];
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
          for (int g = 0; g < len_[2]; ++g) {
            file >> dummy;
          }
        }
      }

  file.close();
}

Real HeliosCK::GetAttenuation(Real g1, Real g2, AirParcel const& var) const {
  // temperature, log-pressure, wave-scaled g-point
  Real val, coord[3] = {var.q[IDN], log(var.q[IPR]), g1};
  interpn(&val, coord, kcoeff_.data(), axis_.data(), len_, 3, 1);
  Real dens = var.q[IPR] / (Constants::kBoltz * var.q[IDN]);
  return exp(val) * dens;  // ln(cm^2 / molecule) -> 1/m
}
