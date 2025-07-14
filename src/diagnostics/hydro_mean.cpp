// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/reconstruct/interpolation.hpp>

// snap
#include <snap/eos/ideal_moist.hpp>

// exchanger
#include <exchanger/exchanger.hpp>

// snap
#include <snap/stride_iterator.hpp>

// canoe
#include <interface/eos.hpp>

// canoe
#include "diagnostics.hpp"

HydroMean::HydroMean(MeshBlock *pmb) : Diagnostics(pmb, "mean"), ncycle_(0) {
  type = "VECTORS";
  std::string varname = "rho_bar,";

  for (int n = 1; n <= NVAPOR; ++n)
    varname += "q" + std::to_string(n) + "_bar,";
  for (int n = 0; n < 3; ++n)
    varname += "vel" + std::to_string(n + 1) + "_bar,";
  varname += "T_bar";

  SetName(varname);
  data.NewAthenaArray(NHYDRO, ncells3_, ncells2_, ncells1_);
}

void HydroMean::Progress(MeshBlock *pmb) {
  Coordinates *pcoord = pmb->pcoord;
  auto const &w = pmb->phydro->w;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  auto peos = pmb->pimpl->peos;

  if (ncycle_ == 0) {
    std::fill(data.data(), data.data() + data.GetSize(), 0.);
  }

  for (int n = 0; n < IPR; ++n)
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          data(n, k, j, i) += w(n, k, j, i);
        }

  // Temperature field is save in IPR
  auto temp = get_temp(peos, w);

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        data(IPR, k, j, i) += temp.at(k, j, i)[0];
      }

  ncycle_++;
}

void HydroMean::Finalize(MeshBlock *pmb) {
  Coordinates *pcoord = pmb->pcoord;
  auto const &w = pmb->phydro->w;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  auto peos = pmb->pimpl->peos;
  auto temp = get_temp(peos, w);

  // hydro mean
  if (ncycle_ > 0) {
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          for (int i = is; i <= ie; ++i) {
            data(n, k, j, i) /= ncycle_;
          }
  } else {
    for (int n = 0; n < IPR; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          for (int i = is; i <= ie; ++i) {
            data(n, k, j, i) = w(n, k, j, i);
          }

    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          data(IPR, k, j, i) = temp.at(k, j, i)[0];
        }
  }

  // clear cycle
  ncycle_ = 0;
}
