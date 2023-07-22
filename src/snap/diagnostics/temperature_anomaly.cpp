#include "../coordinates/coordinates.hpp"
#include "../globals.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "diagnostics.hpp"

TemperatureAnomaly::TemperatureAnomaly(MeshBlock *pmb)
    : Diagnostics(pmb, "tempa") {
  type = "SCALARS";
  long_name = "Temperature anomaly";
  units = "K";
  mean_.NewAthenaArray(ncells1_);
  data.NewAthenaArray(ncells3_, ncells2_, ncells1_);
}

void TemperatureAnomaly::Finalize(AthenaArray<Real> const &w) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pcoord = pmb->pcoord;
  Thermodynamics *pthermo = pmb->pthermo;
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  mean_.ZeroClear();

  // calculate horizontal mean
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pcoord->CellVolume(k, j, is, ie, vol_);
      for (int i = is; i <= ie; ++i)
        mean_(i) += vol_(i) * pthermo->GetTemp(w.at(k, j, i));
    }

  gatherAllData23_(total_vol_, mean_);

  // temperature anomaly
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i)
        data(k, j, i) =
            pthermo->GetTemp(w.at(k, j, i)) - mean_(i) / total_vol_(i);
}
