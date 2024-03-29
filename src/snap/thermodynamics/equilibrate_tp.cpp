// canoe
#include <air_parcel.hpp>

// snap
#include "atm_thermodynamics.hpp"

void Thermodynamics::EquilibrateTP(AirParcel* qfrac) const {
  set_total_equivalent_vapor(qfrac, cloud_index_set_.data(),
                             cloud_reaction_map_);

  // vapor <=> cloud
  for (int i = 1; i <= NVAPOR; ++i) {
    auto rates = TryEquilibriumTP_VaporCloud(*qfrac, i);

    // vapor condensation rate
    qfrac->w[i] += rates[0];

    // cloud concentration rates
    for (int n = 1; n < rates.size(); ++n)
      qfrac->c[cloud_index_set_[i][n - 1]] += rates[n];
  }

  // vapor + vapor <=> cloud
  for (auto const& [ij, info] : cloud_reaction_map_) {
    auto rates = TryEquilibriumTP_VaporVaporCloud(*qfrac, ij);
    auto& indx = info.first;

    // vapor condensation rate
    qfrac->w[indx[0]] += rates[0];
    qfrac->w[indx[1]] += rates[1];

    // cloud concentration rates
    qfrac->c[indx[2]] += rates[3];
  }
}
