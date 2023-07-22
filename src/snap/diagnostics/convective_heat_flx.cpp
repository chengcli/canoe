// C/C++ headers
// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#include "../coordinates/coordinates.hpp"
#include "../reconstruct/interpolation.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "diagnostics.hpp"

ConvectiveHeatFlx::ConvectiveHeatFlx(MeshBlock *pmb)
    : Diagnostics(pmb, "conv_heat_flx") {
  type = "VECTORS";
  grid = "--C";
  long_name = "Z-coordinate convective heat flux (eddy and mean)";
  units = "W/(m^2)";
  eddy_.NewAthenaArray(NHYDRO, ncells3_, ncells2_, ncells1_);
  mean_.NewAthenaArray(NHYDRO, ncells3_, ncells2_, ncells1_);
  // eddy heating rate and mean heating rate
  data.NewAthenaArray(4, 1, 1, ncells1_);
}

ConvectiveHeatFlx::~ConvectiveHeatFlx() {
  eddy_.DeleteAthenaArray();
  mean_.DeleteAthenaArray();
  data.DeleteAthenaArray();
}

void ConvectiveHeatFlx::Progress(AthenaArray<Real> const &w) {
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

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pcoord->CellVolume(k, j, is, ie, vol_);
      for (int i = is; i <= ie; ++i) {
        Real Cp = pthermo->GetMeanCp(w.at(k, j, i));
        data(0, i) += Cp * mean_(IDN, k, j, i) * eddy_(IPR, k, j, i) *
                      eddy_(IVX, k, j, i) * vol_(i);
        data(1, i) += Cp * mean_(IDN, k, j, i) * mean_(IPR, k, j, i) *
                      mean_(IVX, k, j, i) * vol_(i);
        data(2, i) += -mean_(IDN, k, j, i) * mean_(IVX, k, j, i) * 24.79;
        data(3, i) += -eddy_(IDN, k, j, i) * eddy_(IVX, k, j, i) * 24.79;
      }
    }

  ncycle++;
  delete[] data_sum;
}

void ConvectiveHeatFlx::Finalize(AthenaArray<Real> const &w) {
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
    for (int n = 0; n < 4; ++n)
      for (int i = is; i <= ie; ++i) data(n, i) /= ncycle * total_vol[i];
  }

  // clear cycle;
  delete[] total_vol;
  ncycle = 0;
}
