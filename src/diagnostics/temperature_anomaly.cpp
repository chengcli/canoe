// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// diagnostics
#include "diagnostics.hpp"

TemperatureAnomaly::TemperatureAnomaly(MeshBlock *pmb)
    : Diagnostics(pmb, "tempa") {
  type = "SCALARS";
  data.NewAthenaArray(ncells3_, ncells2_, ncells1_);
  mean_.resize(ncells1_);

  pexh = std::make_shared<PlanarExchanger<Real, 2>>("diag/tempa");

  // volumn
  pexh->send_buffer[0].resize(total_vol_.size());
  pexh->recv_buffer[0].resize(total_vol_.size());

  // value
  pexh->send_buffer[1].resize(ncells1_);
  pexh->recv_buffer[1].resize(ncells1_);
}

void TemperatureAnomaly::Finalize(MeshBlock *pmb) {
  Coordinates *pcoord = pmb->pcoord;
  auto const &w = pmb->phydro->w;

  auto pthermo = Thermodynamics::GetInstance();

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  std::fill(total_vol_.begin(), total_vol_.end(), 0.);
  std::fill(mean_.begin(), mean_.end(), 0.);

  // calculate horizontal mean
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pcoord->CellVolume(k, j, is, ie, vol_);

      for (int i = is; i <= ie; ++i) {
        total_vol_[i] += vol_(i);
        mean_[i] += vol_(i) * pthermo->GetTemp(pmb, k, j, i);
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
        data(k, j, i) =
            pthermo->GetTemp(pmb, k, j, i) - mean_[i] / total_vol_[i];
      }
}

void TemperatureAnomaly::packData(MeshBlock const *pmb) {
  pexh->send_buffer[0].swap(total_vol_);
  pexh->send_buffer[1].swap(mean_);
}

void TemperatureAnomaly::unpackData(MeshBlock const *pmb) {
  total_vol_.swap(pexh->recv_buffer[0]);
  mean_.swap(pexh->recv_buffer[1]);
}
