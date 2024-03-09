// C/C++ headers
// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/reconstruct/interpolation.hpp>
#include <athena/thermodynamics/thermodynamics.hpp>

// canoe
#include "diagnostics.hpp"

EddyFlux::EddyFlux(MeshBlock *pmb)
    : Diagnostics(pmb, "eddyflux", "Z-coordinate eddy flux");
{
  type = "VECTORS";
  grid = "--C";
  eddy_.NewAthenaArray(NHYDRO, ncells3_, ncells2_, ncells1_);
  mean_.NewAthenaArray(NHYDRO, ncells3_, ncells2_, ncells1_);
  data.NewAthenaArray(2 * NHYDRO, 1, 1, ncells1_);
}

EddyFlux::~EddyFlux() {
  eddy_.DeleteAthenaArray();
  mean_.DeleteAthenaArray();
  data.DeleteAthenaArray();
}

void EddyFlux::Progress(AthenaArray<Real> const &w) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pcoord = pmb->pcoord;
  Thermodynamics *pthermo = pmb->pthermo;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  Real *data_sum = new Real[(NHYDRO + 1) * ncells1_];
  std::fill(data_sum, data_sum + (NHYDRO + 1) * ncells1_, 0.);

  // calculate horizontal mean
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pcoord->CellVolume(k, j, is, ie, vol_);
      for (int i = is; i <= ie; ++i) {
        for (int n = 0; n < IPR; ++n)
          data_sum[n * ncells1_ + i] += vol_(i) * w(n, k, j, i);
        data_sum[IPR * ncells1_ + i] +=
            vol_(i) * pthermo->GetTemp(w.at(k, j, i));
        data_sum[NHYDRO * ncells1_ + i] += vol_(i);
      }
    }

#ifdef MPI_PARALLEL
  // sum over all ranks
  MPI_Allreduce(MPI_IN_PLACE, data_sum, (NHYDRO + 1) * ncells1_,
                MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif

  // calculate eddy term
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real vol = data_sum[NHYDRO * ncells1_ + i];
        for (int n = 0; n < IPR; ++n) {
          mean_(n, k, j, i) = data_sum[n * ncells1_ + i] / vol;
          eddy_(n, k, j, i) = w(n, k, j, i) - mean_(n, k, j, i);
        }
        mean_(IPR, k, j, i) = data_sum[IPR * ncells1_ + i] / vol;
        eddy_(IPR, k, j, i) =
            pthermo->GetTemp(w.at(k, j, i)) - mean_(IPR, k, j, i);
      }

  // take horizontal average
  if (ncycle == 0) std::fill(data.data(), data.data() + data.GetSize(), 0.);

  for (int n = 0; n < NHYDRO; ++n)
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        pcoord->CellVolume(k, j, is, ie, vol_);
        for (int i = is; i <= ie; ++i) {
          data(n, i) += eddy_(n, k, j, i) * eddy_(IVX, k, j, i) * vol_(i);
          data(n + NHYDRO, i) +=
              mean_(n, k, j, i) * mean_(IVX, k, j, i) * vol_(i);
        }
      }

  ncycle++;
  delete[] data_sum;
}

void EddyFlux::Finalize(AthenaArray<Real> const &w) {
  MeshBlock *pmb = pmy_block_;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  Real *total_vol = new Real[ncells1_];
  std::fill(total_vol, total_vol + ncells1_, 0.);

  // calculate total volume
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pmb->pcoord->CellVolume(k, j, is, ie, vol_);
      for (int i = is; i <= ie; ++i) total_vol[i] += vol_(i);
    }

  // take time and spatial average
  if (ncycle > 0) {
    // sum over all ranks
#ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, data.data(), data.GetSize(), MPI_ATHENA_REAL,
                  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, total_vol, ncells1_, MPI_ATHENA_REAL, MPI_SUM,
                  MPI_COMM_WORLD);
#endif
    for (int n = 0; n < 2 * NHYDRO; ++n)
      for (int i = is; i <= ie; ++i) data(n, i) /= ncycle * total_vol[i];
  }

  // clear cycle;
  delete[] total_vol;
  ncycle = 0;
}
