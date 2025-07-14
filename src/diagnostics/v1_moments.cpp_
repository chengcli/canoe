// athena++
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// exchanger
#include <exchanger/exchanger.hpp>

// canoe
#include "diagnostics.hpp"

V1Moments::V1Moments(MeshBlock *pmb) : Diagnostics(pmb, "w_avg,w2_avg,w3_avg") {
  type = "VECTORS";

  data.NewAthenaArray(3, ncells1_);

  pexh = std::make_shared<PlanarExchanger<Real, 2>>("diag/v1_moments");

  // volumn
  pexh->send_buffer[0].resize(total_vol_.size());
  pexh->recv_buffer[0].resize(total_vol_.size());
  // value
  pexh->send_buffer[1].resize(3 * ncells1_);
  pexh->recv_buffer[1].resize(3 * ncells1_);
}

void V1Moments::Finalize(MeshBlock *pmb) {
  std::fill(total_vol_.begin(), total_vol_.end(), 0.);
  std::fill(data.data(), data.data() + data.GetSize(), 0.);

  // calculate total volume
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol_);

      for (int i = pmb->is; i <= pmb->ie; ++i) {
        total_vol_[i] += vol_(i);
        Real v1 = pmb->phydro->w(IVX, k, j, i);
        data(0, i) += vol_(i) * v1;
        data(1, i) += vol_(i) * v1 * v1;
        data(2, i) += vol_(i) * v1 * v1 * v1;
      }
    }

  pexh->Regroup(pmb, X1DIR);
  packData(pmb);
  pexh->ReduceAll(pmb, MPI_SUM);
  unpackData(pmb);

  for (int i = pmb->is; i <= pmb->ie; ++i)
    for (int n = 0; n < 3; ++n) {
      data(n, i) /= total_vol_[i];
    }
}

void V1Moments::packData(MeshBlock const *pmb) {
  pexh->send_buffer[0].swap(total_vol_);

  for (int n = 0; n < 3; ++n)
    for (int i = pmb->is; i <= pmb->ie; ++i) {
      pexh->send_buffer[1][n * ncells1_ + i] = data(n, i);
    }
}

void V1Moments::unpackData(MeshBlock const *pmb) {
  total_vol_.swap(pexh->recv_buffer[0]);

  for (int n = 0; n < 3; ++n) {
    for (int i = pmb->is; i <= pmb->ie; ++i) {
      data(n, i) = pexh->recv_buffer[1][n * ncells1_ + i];
    }
  }
}
