// C/C++
#include <sstream>

// athena
#include <athena/athena.hpp>

// application
#include <application/exceptions.hpp>

// canoe
#include <variable.hpp>

// climath
#include <climath/core.h>

#include <climath/root.hpp>

// snap
#include "thermodynamics.hpp"

// Calculates phase equilibrium of
// a * Vapor1 + b * Vapor2 = c * Cloud
//
// Example phase equilibrium:
// NH3 + H2S -> NH4SH
//
RealArray3 Thermodynamics::TryEquilibriumTP_VaporVaporCloud(Variable const& air,
                                                            IndexPair ij,
                                                            Real cv_hat,
                                                            bool misty) const {
  auto& info = cloud_reaction_map_.at(ij);

  int j1 = info.first[0];
  int j2 = info.first[1];
  int j3 = info.first[2];
  auto& stoi = info.second;

  Real a = stoi[0] / stoi[2], b = stoi[1] / stoi[2];

  Real ptol = air.w[IPR], temp = air.w[IDN], gmol = 1.;

#pragma omp simd reduction(+ : gmol)
  for (int n = 0; n < NCLOUD; ++n) gmol += -air.c[n];

  gmol -= air.w[j1] + air.w[j2];

  Real x0 = a * air.c[j3] + air.w[j1], y0 = b * air.c[j3] + air.w[j2], z0 = 0.,
       pp1 = x0 / (x0 + y0 + gmol) * ptol, pp2 = y0 / (x0 + y0 + gmol) * ptol;

  Real mol0 = air.w[j1], mol1 = air.w[j2];
  Real svp = svp_func2_.at(ij)(air, j1, j2, j3);

  if (pp1 * pp2 > svp) {
    Real zeta = svp / (ptol * ptol), rh;

    int error;

    if (a * y0 > b * x0) {
      error = root(0., 1., 1.E-8, &rh, [a, b, x0, y0, gmol, zeta](double r) {
        Real x = r * x0, y = y0 - b / a * (x0 - x);
        return (x * y) / sqr(x + y + gmol) - zeta;
      });
      y0 -= b / a * (1. - rh) * x0;
      z0 = (1. - rh) * x0 / a;
      x0 *= rh;
    } else {
      error = root(0., 1., 1.E-8, &rh, [a, b, x0, y0, gmol, zeta](double r) {
        Real y = r * y0, x = x0 - a / b * (y0 - y);
        return (x * y) / sqr(x + y + gmol) - zeta;
      });
      x0 -= a / b * (1. - rh) * y0;
      z0 = (1. - rh) * y0 / b;
      y0 *= rh;
    }

    // correction for y0 << x0 and x0 >> y0
    if (y0 < 1.E-8 * x0) {
      y0 = (x0 + gmol) * (x0 + gmol) / x0 * zeta;
      z0 = (b * air.c[j3] + air.w[j2] - y0) / b;
      x0 = a * air.c[j3] + air.w[j1] - a * z0;
    } else if (x0 < 1.E-8 * y0) {
      x0 = (y0 + gmol) * (y0 + gmol) / y0 * zeta;
      z0 = (a * air.c[j3] + air.w[j1] - x0) / a;
      y0 = b * air.c[j3] + air.w[j2] - b * z0;
    }

    if (error) {
      std::stringstream msg;
      msg << "Condensation a(A) + b(B) -> c(C) failed" << std::endl;
      msg << "Air parcel = (" << air << ")" << std::endl;
      throw RuntimeError("TryEquilibriumTP_VaporVaporCloud", msg.str());
    }
  }

  std::array<Real, 3> rates;
  rates[0] = x0 - air.w[j1];
  rates[1] = y0 - air.w[j2];
  rates[2] = z0 - air.c[j3];

  return rates;
}
