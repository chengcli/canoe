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
#include <air_parcel.hpp>
#include <constants.hpp>

// snap
#include "atm_thermodynamics.hpp"

void Thermodynamics::SaturationAdjustment(AirColumn& air_column) const {
  // return if there's no vapor
  if (NVAPOR == 0) return;

  for (auto& air : air_column) {
    air.ToMoleFraction();
    AirParcel air0(air);

    // calculate total mole density
    Real rmole = get_density_mole(air);
    Real umole = get_internal_energy_mole(air, cv_ratio_mole_.data(),
                                          latent_energy_mole_.data(),
                                          cloud_index_set_.data());

    int iter = 0;

    Real Teq = air.w[IDN];
    std::stringstream msg;

    while (iter++ < sa_max_iter_) {
      // msg << "iter = " << iter << std::endl;
      // msg << "old var = " << air << std::endl;

      Real cvd = 1. / (get_gammad(air) - 1.);
      Real fsig = 1.;

      bool isNeg = false;

      for (int i = 1; i <= NVAPOR; ++i) {
        Real qv = air.w[i];
        if (qv < 0) {
          isNeg = true;
          std::cout << "gas " << i << " is negative: " << qv << std::endl;
          std::cout << "Temp is " << air.w[IDN] << std::endl;
        }

#pragma omp simd reduction(+ : qv)
        for (auto j : cloud_index_set_[i]) qv += air.c[j];

        fsig += (cv_ratio_mole_[i] - 1.) * qv;
      }

      // condensation
      for (int i = 1; i <= NVAPOR; ++i) {
        auto rates = TryEquilibriumTP_VaporCloud(air, i, cvd * fsig);

        // vapor condensation rate
        air.w[i] += rates[0];
        if (isNeg && i == 1) {
          std::cout << "new qv: " << air.w[i] << std::endl;
          std::cout << "rate: " << rates[0] << std::endl;
        }
        // cloud concentration rates
        for (int j = 1; j < rates.size(); ++j) {
          air.c[cloud_index_set_[i][j - 1]] += rates[j];
        }
      }

      Real Told = air.w[IDN];
      update_TP_conserving_U(&air, rmole, umole, cv_ratio_mole_.data(),
                             latent_energy_mole_.data(),
                             cloud_index_set_.data());
      // msg << "new var = " << air << std::endl;
      if (fabs(air.w[IDN] - Teq) < sa_ftol_) break;

      // relax temperature and pressure
      Teq = air.w[IDN];
      Real qgas = 1.;
#pragma omp simd reduction(+ : qgas)
      for (int n = 0; n < NCLOUD; ++n) qgas += -air.c[n];

      air.w[IDN] = Told + sa_relax_ * (air.w[IDN] - Told);
      air.w[IPR] = rmole * qgas * Constants::Rgas * air.w[IDN];
    }

    if (iter > sa_max_iter_) {
      msg << "Variables before iteration q0 = (" << air0 << ")" << std::endl;
      msg << "Variables after iteration q = (" << air << ")" << std::endl;
      throw RuntimeError("SaturationAdjustment", msg.str());
    }
  }
}
