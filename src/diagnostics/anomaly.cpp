// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>

// snap
#include <snap/eos/ideal_moist.hpp>

// snap
#include <snap/stride_iterator.hpp>

// exchanger
#include <exchanger/exchanger.hpp>

// canoe
#include <interface/eos.hpp>

// diagnostics
#include "diagnostics.hpp"

Anomaly::Anomaly(MeshBlock *pmb) : Diagnostics(pmb, "rhoa,tempa,v1a,presa") {
  type = "VECTORS";
  data.NewAthenaArray(4, ncells3_, ncells2_, ncells1_);
  mean_.resize(4 * ncells1_);

  pexh = std::make_shared<PlanarExchanger<Real, 2>>("diag/anomaly");

  // volumn
  pexh->send_buffer[0].resize(total_vol_.size());
  pexh->recv_buffer[0].resize(total_vol_.size());

  // value
  pexh->send_buffer[1].resize(4 * ncells1_);
  pexh->recv_buffer[1].resize(4 * ncells1_);
}

void Anomaly::Finalize(MeshBlock *pmb) {
  Coordinates *pcoord = pmb->pcoord;
  auto const &w = pmb->phydro->w;

  auto peos = pmb->pimpl->peos;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  std::fill(total_vol_.begin(), total_vol_.end(), 0.);
  std::fill(mean_.begin(), mean_.end(), 0.);

  // calculate horizontal mean
  auto temp = get_temp(peos, w);

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pcoord->CellVolume(k, j, is, ie, vol_);

      for (int i = is; i <= ie; ++i) {
        total_vol_[i] += vol_(i);
        // density anomaly
        mean_[i] += vol_(i) * w(IDN, k, j, i);
        // temperature anomaly
        mean_[ncells1_ + i] += vol_(i) * temp.at(k, j, i)[0];
        // vertical velocity anomaly
        mean_[2 * ncells1_ + i] += vol_(i) * w(IVX, k, j, i);
        // pressure anomaly
        mean_[3 * ncells1_ + i] += vol_(i) * w(IPR, k, j, i);
      }
    }

  pexh->Regroup(pmb, X1DIR);
  packData(pmb);
  pexh->ReduceAll(pmb, MPI_SUM);
  unpackData(pmb);

  // pressure anomaly
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        // density anomaly
        data(0, k, j, i) = w(IDN, k, j, i) - mean_[i] / total_vol_[i];
        // temperature anomaly
        data(1, k, j, i) =
            temp.at(k, j, i)[0] - mean_[ncells1_ + i] / total_vol_[i];
        // vertical velocity anomaly
        data(2, k, j, i) =
            w(IVX, k, j, i) - mean_[2 * ncells1_ + i] / total_vol_[i];
        // pressure anomaly
        data(3, k, j, i) =
            w(IPR, k, j, i) - mean_[3 * ncells1_ + i] / total_vol_[i];
      }
}

void Anomaly::packData(MeshBlock const *pmb) {
  pexh->send_buffer[0].swap(total_vol_);
  pexh->send_buffer[1].swap(mean_);
}

void Anomaly::unpackData(MeshBlock const *pmb) {
  total_vol_.swap(pexh->recv_buffer[0]);
  mean_.swap(pexh->recv_buffer[1]);
}
