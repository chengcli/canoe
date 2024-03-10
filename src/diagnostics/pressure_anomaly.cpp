// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>

// canoe
#include "diagnostics.hpp"

PressureAnomaly::PressureAnomaly(MeshBlock *pmb) : Diagnostics(pmb, "presa") {
  type = "SCALARS";
  mean_.NewAthenaArray(ncells1_);
  data.NewAthenaArray(ncells3_, ncells2_, ncells1_);

  Regroup(pmb, X1DIR);

  // calculate total volume
  total_vol_.ZeroClear();
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol_);

      for (int i = pmb->is; i <= pmb->ie; ++i) {
        total_vol_(i) += vol_(i);
      }
    }

#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, total_vol_.data(), total_vol_.GetSize(),
                MPI_ATHENA_REAL, MPI_SUM, mpi_comm_);
#endif
}

void PressureAnomaly::Finalize(MeshBlock *pmb) {
  Coordinates *pcoord = pmb->pcoord;
  auto const &w = pmb->phydro->w;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  mean_.ZeroClear();

  // calculate horizontal mean
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pcoord->CellVolume(k, j, is, ie, vol_);

      for (int i = is; i <= ie; ++i) {
        mean_(i) += vol_(i) * w(IPR, k, j, i);
      }
    }

    // gatherAllData23_(total_vol_, mean_);
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, mean_.data(), mean_.GetSize(), MPI_ATHENA_REAL,
                MPI_SUM, mpi_comm_);
#endif

  // pressure anomaly
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        data(k, j, i) = w(IPR, k, j, i) - mean_(i) / total_vol_(i);
      }
}
