// C/C++ headers
// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#include "../coordinates/coordinates.hpp"
#include "../radiation/radiation.hpp"
#include "diagnostics.hpp"

RadiativeFlux::RadiativeFlux(MeshBlock *pmb) : Diagnostics(pmb, "radflux") {
  type = "VECTORS";
  grid = "--F";
  long_name = "total upward radiative flux,total downward radiative flux";
  units = "w/m^2";
  // 0: upward flux
  // 1: downward flux
  data.NewAthenaArray(2, 1, 1, ncells1_ + 1);
}

void RadiativeFlux::Progress(AthenaArray<Real> const &w) {
  MeshBlock *pmb = pmy_block_;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  if (ncycle == 0) std::fill(data.data(), data.data() + data.GetSize(), 0.);

  // sum over horizontal grids weighted by area
  RadiationBand *pband = pmb->prad->pband;
  while (pband != NULL) {
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        pmb->pcoord->Face1Area(k, j, is, ie + 1, x1area_);
        for (int i = is; i <= ie + 1; ++i) {
          data(0, i) += x1area_(i) * pband->bflxup(k, j, i);
          data(1, i) += x1area_(i) * pband->bflxdn(k, j, i);
        }
      }
    pband = pband->next;
  }

  ncycle++;
}

void RadiativeFlux::Finalize(AthenaArray<Real> const &w) {
  MeshBlock *pmb = pmy_block_;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  Real *total_area = new Real[ncells1_ + 1];
  std::fill(total_area, total_area + ncells1_ + 1, 0.);

  // calculate total face area
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pmb->pcoord->Face1Area(k, j, is, ie + 1, x1area_);
      for (int i = is; i <= ie + 1; ++i) total_area[i] += x1area_(i);
    }

  // take time and spatial average
  if (ncycle > 0) {
    // sum over all ranks
#ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, data.data(), data.GetSize(), MPI_ATHENA_REAL,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, total_area, ncells1_ + 1, MPI_ATHENA_REAL,
                  MPI_SUM, MPI_COMM_WORLD);
#endif
    for (int i = is; i <= ie + 1; ++i) {
      data(0, i) /= ncycle * total_area[i];
      data(1, i) /= ncycle * total_area[i];
    }
  }

  // clear cycle;
  delete[] total_area;
  ncycle = 0;
}
