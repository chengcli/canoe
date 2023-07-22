// C/C++
#include <cstdlib>
#include <cstring>
#include <sstream>

// athena
#include <athena/athena_arrays.hpp>
#include <athena/hydro/hydro.hpp>

// application
#include <application/exceptions.hpp>

// canoe
#include <variable.hpp>

// snap
#include "thermodynamics.hpp"

void Thermodynamics::SaturationAdjustment(Variable *qfrac) const {
  // return if there's no vapor
  if (NVAPOR == 0) return;
  Variable qfrac0(*qfrac);

  // calculate total mole density
  Real rmole = getDensityMole(*qfrac);
  Real umole = getInternalEnergyMole(*qfrac);

  int iter = 0;

  Real Teq = qfrac->w[IDN];
  while (iter++ < sa_max_iter_) {
    Real cvd = Rd_ / (GetGammad(*qfrac) - 1.);
    Real fsig = 1.;
    for (int i = 1; i <= NVAPOR; ++i) {
      Real qv = qfrac->w[i];

#pragma omp simd reduction(+ : qv)
      for (auto j : cloud_index_set_[i]) qv += qfrac->c[j];

      fsig += (cv_ratio_mole_[i] - 1.) * qv;
    }

    // condensation
    for (int i = 1; i <= NVAPOR; ++i) {
      auto rates = TryEquilibriumTP_VaporCloud(*qfrac, i, cvd * fsig);

      // vapor condensation rate
      qfrac->w[i] += rates[0];

      // cloud concentration rates
      for (int j = 1; j < rates.size(); ++j) {
        qfrac->c[cloud_index_set_[i][j - 1]] += rates[j];
      }
    }

    Real Told = qfrac->w[IDN];
    updateTPConservingU(qfrac, rmole, umole);
    if (fabs(qfrac->w[IDN] - Teq) < sa_ftol_) break;

    // relax temperature and pressure
    Teq = qfrac->w[IDN];
    Real qgas = 1.;
#pragma omp simd reduction(+ : qgas)
    for (int n = 0; n < NCLOUD; ++n) qgas += -qfrac->c[n];

    qfrac->w[IDN] = Told + sa_relax_ * (qfrac->w[IDN] - Told);
    qfrac->w[IPR] = rmole * qgas * Constants::Rgas * qfrac->w[IDN];
  }

  if (iter > sa_max_iter_) {
    std::stringstream msg;
    msg << "Variables before iteration q0 = (" << qfrac0 << ")" << std::endl;
    msg << "Variables after iteration q = (" << *qfrac << ")" << std::endl;
    throw RuntimeError("SaturationAdjustment", msg.str());
  }
}
