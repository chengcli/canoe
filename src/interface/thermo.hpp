#pragma once

// athena
#include <athena/athena.hpp>

// kintera
#include <kintera/constants.h>

#include <kintera/thermo/thermo.hpp>

inline static Real get_rd() {
  return kintera::constants::Rgas / kintera::species_weights[0];
}

inline static Real get_gammad() {
  auto cvd = kintera::species_cref_R[0];
  return (cvd + 1.) / cvd;
}

inline Real get_cv_ratio(kintera::ThermoY const& pthermo, int n) {
  int nvapor = pthermo->options.vapor_ids().size();
  int ncloud = pthermo->options.cloud_ids().size();

  if (n < 0 || n >= nvapor + ncloud) {
    throw std::out_of_range("Invalid species index: " + std::to_string(n));
  }

  int id;
  if (n < nvapor - 1) {  // vapor
    id = pthermo->options.vapor_ids()[n + 1];
  } else {  // cloud
    id = pthermo->options.cloud_ids()[n - nvapor + 1];
  }

  auto cv_ratio = kintera::species_cref_R[id] / kintera::species_cref_R[0];
  auto mu_ratio = kintera::species_weights[id] / kintera::species_weights[0];
  return cv_ratio / mu_ratio;
}

inline Real get_inv_mu_ratio(kintera::ThermoY const& pthermo, int n) {
  int nvapor = pthermo->options.vapor_ids().size();
  int ncloud = pthermo->options.cloud_ids().size();

  if (n < 0 || n >= nvapor + ncloud) {
    throw std::out_of_range("Invalid species index: " + std::to_string(n));
  }

  int id;
  if (n < nvapor - 1) {  // vapor
    id = pthermo->options.vapor_ids()[n + 1];
  } else {  // cloud
    id = pthermo->options.cloud_ids()[n - nvapor + 1];
  }

  return kintera::species_weights[0] / kintera::species_weights[id];
}
