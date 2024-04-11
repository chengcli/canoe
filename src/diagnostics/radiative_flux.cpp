// athena
#include <athena/coordinates/coordinates.hpp>

// exchanger
#include <exchanger/exchanger.hpp>

// harp
#include <harp/radiation.hpp>

// diagnostics
#include "diagnostics.hpp"

RadiativeFlux::RadiativeFlux(MeshBlock *pmb)
    : Diagnostics(pmb, "rflx_up,rflx_dn"), ncycle_(0) {
  type = "VECTORS";

  // 0: upward flux
  // 1: downward flux
  data.NewAthenaArray(2, ncells1_ + 1);

  pexh = std::make_shared<PlanarExchanger<Real, 2>>("diag/rflx");

  // volumn
  pexh->send_buffer[0].resize(total_area_.size());
  pexh->recv_buffer[0].resize(total_area_.size());

  // value
  pexh->send_buffer[1].resize(2 * (ncells1_ + 1));
  pexh->recv_buffer[1].resize(2 * (ncells1_ + 1));
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

    // calculate total face area
    std::fill(total_area_.begin(), total_area_.end(), 0.);
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j) {
        pmb->pcoord->Face1Area(k, j, pmb->is, pmb->ie + 1, x1area_);

        for (int i = pmb->is; i <= pmb->ie + 1; ++i) {
          total_area_[i] += x1area_(i);
        }
      }

    pexh->Regroup(pmb, X1DIR);

    packData(pmb);
    pexh->ReduceAll(pmb, MPI_SUM);
    unpackData(pmb);

    for (int i = pmb->is; i <= pmb->ie + 1; ++i) {
      data(0, i) /= ncycle_ * total_area_[i];
      data(1, i) /= ncycle_ * total_area_[i];
    }
  }

  // clear cycle;
  ncycle_ = 0;
}

void RadiativeFlux::packData(MeshBlock const *pmb) {
  pexh->send_buffer[0].swap(total_area_);

  for (int n = 0; n < 2; ++n) {
    for (int i = pmb->is; i <= pmb->ie + 1; ++i) {
      pexh->send_buffer[1][n * (ncells1_ + 1) + i] = data(n, i);
    }
  }
}

void RadiativeFlux::unpackData(MeshBlock const *pmb) {
  total_area_.swap(pexh->recv_buffer[0]);

  for (int n = 0; n < 2; ++n) {
    for (int i = pmb->is; i <= pmb->ie + 1; ++i) {
      data(n, i) = pexh->recv_buffer[1][n * (ncells1_ + 1) + i];
    }
  }
}
