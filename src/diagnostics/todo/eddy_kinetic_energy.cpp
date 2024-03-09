/** @file eddy_kinetic_energy.cpp
 * @brief Implement eddy kinetic energy and mean kinetic energy output (eke)
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Tuesday May 18, 2021 08:07:22 PDT
 * @bug No known bugs.
 */

// Athena++ headers
#include <athena/coordinates/coordinates.hpp>
#include <athena/globals.hpp>
#include <athena/utils/utils.hpp>

// canoe
#include "diagnostics.hpp"

// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

EddyKineticEnergy::EddyKineticEnergy(MeshBlock *pmb)
    : Diagnostics(pmb, "eke", "eddy kinetic energy, mean kinetic energy") {
  type = "VECTORS";
  grid = "--C";
  units = "J/m^3,J/m^3";
  data.NewAthenaArray(2, 1, 1, ncells1_);
}

void EddyKineticEnergy::Finalize(AthenaArray<Real> const &w) {
  MeshBlock *pmb = pmy_block_;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  Real *data_sum = new Real[4 * ncells1_];
  Real **dir_sum;
  NewCArray(dir_sum, 3, 4 * ncells1_);

  std::fill(data_sum, data_sum + 4 * ncells1_, 0.);
  // calculate horizontal mean
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pmb->pcoord->CellVolume(k, j, is, ie, vol_);
      for (int i = is; i <= ie; ++i) {
        data_sum[0 * ncells1_ + i] += vol_(i) * w(IVX, k, j, i);
        data_sum[1 * ncells1_ + i] += vol_(i) * w(IVY, k, j, i);
        data_sum[2 * ncells1_ + i] += vol_(i) * w(IVZ, k, j, i);
        data_sum[3 * ncells1_ + i] += vol_(i);
      }
    }

  /*#ifdef MPI_PARALLEL
    MPI_Comm comm;

          for (int d = 0; d < 3; ++d) {
                  if ((d == 1) && (!pmb->pmy_mesh->f2)) break;
                  if ((d == 2) && (!pmb->pmy_mesh->f3)) break;
                  SetColor(static_cast<CoordinateDirection>(d));
                  if (Globals::my_rank == 0) {
                          for (int i = 0; i < Globals::nranks; ++i)
                                  std::cout << color_[i] << " ";
                          std::cout << std::endl;
                  }
                  MPI_Comm_split(MPI_COMM_WORLD, color_[Globals::my_rank],
  Globals::my_rank, &comm); MPI_Allreduce(data_sum, dir_sum[d], 4*ncells1_,
  MPI_ATHENA_REAL, MPI_SUM, comm); MPI_Comm_free(&comm);
          }
  #endif*/

  // calculate eddy and mean
  data.ZeroClear();
  for (int n = 0; n < 3; ++n) {
    if ((n == 1) && (!pmb->pmy_mesh->f2)) break;
    if ((n == 2) && (!pmb->pmy_mesh->f3)) break;
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          Real mean =
              dir_sum[n][n * ncells1_ + i] / dir_sum[n][3 * ncells1_ + i];
          Real eddy = w(IVX + n, k, j, i) - mean;
          data(0, i) += 0.5 * w(IDN, k, j, i) * eddy * eddy;
          data(1, i) += 0.5 * w(IDN, k, j, i) * mean * mean;
        }
  }

  delete[] data_sum;
  FreeCArray(dir_sum);
}
