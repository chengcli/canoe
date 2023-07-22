#include "../coordinates/coordinates.hpp"
#include "../globals.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "diagnostics.hpp"

HydroFlux::HydroFlux(MeshBlock *pmb) : Diagnostics(pmb, "hydroflux") {
  type = "VECTORS";
  grid = "--C";
  varname = "v1rho,";
  long_name = "horizontally averaged mass flux,";
  units = "kg/(m^2.s),";

  for (int n = 1; n <= NVAPOR; ++n) {
    varname += "v1q" + std::to_string(n) + ",";
    units += "kg/(kg.m^2.s),";
    long_name += "horizontally averaged vapor mass flux,";
  }
  for (int n = 0; n < 3; ++n) {
    units += "m/s^2,";
    varname += "v1v" + std::to_string(n + 1) + ",";
    long_name += "horizontally averaged momentum flux,";
  }
  units += "m.K/s";
  varname += "v1T";
  long_name += "horizontally averaged heat flux";

  data.NewAthenaArray(NHYDRO, 1, 1, ncells1_);
}

void HydroFlux::Progress(AthenaArray<Real> const &w) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pcoord = pmb->pcoord;
  Thermodynamics *pthermo = pmb->pthermo;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  if (ncycle == 0) std::fill(data.data(), data.data() + data.GetSize(), 0.);

  // sum over horizontal grids weighted by volume
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pcoord->CellVolume(k, j, is, ie, vol_);
      for (int i = is; i <= ie; ++i) {
        data(0, i) += vol_(i) * w(IDN, k, j, i) * w(IVX, k, j, i);
        for (int n = 1; n < IPR; ++n)
          data(n, i) += vol_(i) * w(n, k, j, i) * w(IVX, k, j, i);
        data(IPR, i) +=
            vol_(i) * pthermo->GetTemp(w.at(k, j, i)) * w(IVX, k, j, i);
      }
    }

  ncycle++;
}

void HydroFlux::Finalize(AthenaArray<Real> const &w) {
  MeshBlock *pmb = pmy_block_;

  gatherAllData23_(total_vol_, data);

  // take time and spatial average
  if (ncycle > 0) {
    for (int n = 0; n < NHYDRO; ++n)
      for (int i = pmb->is; i <= pmb->ie; ++i)
        data(n, i) /= ncycle * total_vol_(i);
  }

  // clear cycle;
  ncycle = 0;
}
