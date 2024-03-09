// athena
#include <athena/coordinates/coordinates.hpp>

// canoe
#include "diagnostics.hpp"

PressureAnomaly::PressureAnomaly(MeshBlock *pmb)
    : Diagnostics(pmb, "presa", "Pressure anomaly") {
  type = "SCALARS";
  units = "pa";
  mean_.NewAthenaArray(ncells1_);
  data.NewAthenaArray(ncells3_, ncells2_, ncells1_);
}

void PressureAnomaly::Finalize(AthenaArray<Real> const &w) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pcoord = pmb->pcoord;
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  mean_.ZeroClear();

  // calculate horizontal mean
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pcoord->CellVolume(k, j, is, ie, vol_);
      for (int i = is; i <= ie; ++i) mean_(i) += vol_(i) * w(IPR, k, j, i);
    }

  gatherAllData23_(total_vol_, mean_);

  // pressure anomaly
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i)
        data(k, j, i) = w(IPR, k, j, i) - mean_(i) / total_vol_(i);
}
