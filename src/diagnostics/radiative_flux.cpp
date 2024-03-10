// C/C++ headers
// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// athena
#include <athena/coordinates/coordinates.hpp>

// harp
#include <harp/radiation.hpp>

// canoe
#include <configure.hpp>

// diagnostics
#include "diagnostics.hpp"

RadiativeFlux::RadiativeFlux(MeshBlock *pmb)
    : Diagnostics(pmb, "rflx_up,rflx_dn"), ncycle_(0) {
  type = "VECTORS";

  // 0: upward flux
  // 1: downward flux
  data.NewAthenaArray(2, ncells1_ + 1);

  Regroup(pmb, X1DIR);

  // calculate total face area
  total_area_.ZeroClear();
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      pmb->pcoord->Face1Area(k, j, pmb->is, pmb->ie + 1, x1area_);

      for (int i = pmb->is; i <= pmb->ie + 1; ++i) {
        total_area_(i) += x1area_(i);
      }
    }

#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, total_area_.data(), total_area_.GetSize(),
                MPI_ATHENA_REAL, MPI_SUM, mpi_comm_);
#endif
}

void RadiativeFlux::Progress(MeshBlock *pmb) {
  Coordinates *pcoord = pmb->pcoord;
  auto prad = pmb->pimpl->prad;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  if (ncycle_ == 0) {
    std::fill(data.data(), data.data() + data.GetSize(), 0.);
  }

  // sum over horizontal grids weighted by area
  for (int b = 0; b < prad->GetNumBands(); ++b) {
    auto pband = prad->GetBand(b);
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        pcoord->Face1Area(k, j, is, ie + 1, x1area_);
        for (int i = is; i <= ie + 1; ++i) {
          data(0, i) += x1area_(i) * pband->bflxup(k, j, i);
          data(1, i) += x1area_(i) * pband->bflxdn(k, j, i);
        }
      }
  }
  ncycle_++;
}

void RadiativeFlux::Finalize(MeshBlock *pmb) {
  // take time and spatial average
  if (ncycle_ > 0) {
    // sum over all ranks
#ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, data.data(), data.GetSize(), MPI_ATHENA_REAL,
                  MPI_SUM, mpi_comm_);
#endif
    for (int i = pmb->is; i <= pmb->ie + 1; ++i) {
      data(0, i) /= ncycle_ * total_area_(i);
      data(1, i) /= ncycle_ * total_area_(i);
    }
  }

  // clear cycle;
  ncycle_ = 0;
}
