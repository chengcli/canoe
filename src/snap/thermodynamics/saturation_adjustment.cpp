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
#include <snap/thermodynamics/thermodynamics.hpp>

void Thermodynamics::SaturationAdjustment(Variable *qfrac) const {
  // return if there's no vapor
  if (NVAPOR == 0) return;
  Variable qfrac0(*qfrac);

  // calculate total mole density
  Real rmole = getDensityMole(*qfrac);
  Real umole = getInternalEnergyMole(*qfrac);

  int iter = 0;
  while (iter++ < sa_max_iter_) {
    // condensation
    for (int i = 1; i <= NVAPOR; ++i) {
      auto rates = TryEquilibriumTP(*qfrac, i, true);

      // vapor condensation rate
      qfrac->w[i] += rates[0];

      // cloud concentration rates
      for (int j = 1; j < rates.size(); ++j) {
        qfrac->c[cloud_index_set_[i][j - 1]] += rates[j];
      }
    }

    Real Told = qfrac->w[IDN];
    updateTPConservingU(qfrac, rmole, umole);
    if (fabs(qfrac->w[IDN] - Told) < sa_ftol_) break;
  }

  if (iter >= sa_max_iter_) {
    std::stringstream msg;
    msg << "Variables before iteration q0 = (" << qfrac0 << ")" << std::endl;
    msg << "Variables after iteration q = (" << *qfrac << ")" << std::endl;
    throw RuntimeError("SaturationAdjustment", msg.str());
  }
}
