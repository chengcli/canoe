// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>

// snap
#include <snap/eos/ideal_moist.hpp>

// exchanger
#include <exchanger/exchanger.hpp>

// snap
#include <snap/stride_iterator.hpp>

// canoe
#include <interface/eos.hpp>

// diagnostics
#include "diagnostics.hpp"

HydroFlux::HydroFlux(MeshBlock *pmb)
    : Diagnostics(pmb, "hydroflux"), ncycle_(0) {
  type = "VECTORS";
  std::string varname = "v1rho,";

  for (int n = 1; n <= NVAPOR; ++n) varname += "v1q" + std::to_string(n) + ",";
  for (int n = 0; n < 3; ++n) varname += "v1v" + std::to_string(n + 1) + ",";
  varname += "v1T";

  SetName(varname);
  data.NewAthenaArray(NHYDRO, ncells1_);

  pexh = std::make_shared<PlanarExchanger<Real, 2>>("diag/hydroflux");

  // volumn
  pexh->send_buffer[0].resize(total_vol_.size());
  pexh->recv_buffer[0].resize(total_vol_.size());
  // value
  pexh->send_buffer[1].resize(NHYDRO * ncells1_);
  pexh->recv_buffer[1].resize(NHYDRO * ncells1_);
}

void HydroFlux::Progress(MeshBlock *pmb) {
  Coordinates *pcoord = pmb->pcoord;
  auto const &w = pmb->phydro->w;

  auto peos = pmb->pimpl->peos;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  if (ncycle_ == 0) {
    std::fill(data.data(), data.data() + data.GetSize(), 0.);
  }

  auto temp = get_temp(peos, w);

  // sum over horizontal grids weighted by volume
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pcoord->CellVolume(k, j, is, ie, vol_);

      for (int i = is; i <= ie; ++i) {
        data(0, i) += vol_(i) * w(IDN, k, j, i) * w(IVX, k, j, i);
        for (int n = 1; n < IPR; ++n)
          data(n, i) += vol_(i) * w(n, k, j, i) * w(IVX, k, j, i);
        data(IPR, i) += vol_(i) * temp.at(k, j, i)[0] * w(IVX, k, j, i);
      }
    }

  ncycle_++;
}

void HydroFlux::Finalize(MeshBlock *pmb) {
  // take time and spatial average
  if (ncycle_ > 0) {
    std::fill(total_vol_.begin(), total_vol_.end(), 0.);

    // calculate total volume
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j) {
        pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol_);

        for (int i = pmb->is; i <= pmb->ie + 1; ++i) {
          total_vol_[i] += vol_(i);
        }
      }

    pexh->Regroup(pmb, X1DIR);

    packData(pmb);
    pexh->ReduceAll(pmb, MPI_SUM);
    unpackData(pmb);

    for (int n = 0; n < NHYDRO; ++n) {
      for (int i = pmb->is; i <= pmb->ie; ++i)
        data(n, i) /= ncycle_ * total_vol_[i];
    }
  }

  // clear cycle;
  ncycle_ = 0;
}

void HydroFlux::packData(MeshBlock const *pmb) {
  pexh->send_buffer[0].swap(total_vol_);

  for (int n = 0; n < NHYDRO; ++n) {
    for (int i = pmb->is; i <= pmb->ie; ++i) {
      pexh->send_buffer[1][n * ncells1_ + i] = data(n, i);
    }
  }
}

void HydroFlux::unpackData(MeshBlock const *pmb) {
  total_vol_.swap(pexh->recv_buffer[0]);

  for (int n = 0; n < NHYDRO; ++n) {
    for (int i = pmb->is; i <= pmb->ie; ++i) {
      data(n, i) = pexh->recv_buffer[1][n * ncells1_ + i];
    }
  }
}
