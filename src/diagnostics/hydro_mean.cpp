// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/reconstruct/interpolation.hpp>
#include <athena/stride_iterator.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// canoe
#include "diagnostics.hpp"

HydroMean::HydroMean(MeshBlock *pmb) : Diagnostics(pmb, "mean") {
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

HydroMean::~HydroMean() { data.DeleteAthenaArray(); }

void HydroMean::Progress(AthenaArray<Real> const &w) {
  MeshBlock *pmb = pmy_block_;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  auto pthermo = Thermodynamics::GetInstance();

  if (ncycle_ == 0) {
    std::fill(data.data(), data.data() + data.GetSize(), 0.);
  }

  for (int n = 0; n < IPR; ++n)
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) data(n, k, j, i) += w(n, k, j, i);

  // Temperature field is save in IPR
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i)
        data(IPR, k, j, i) += pthermo->GetTemp(pmb, k, j, i);

  ncycle_++;
}

void HydroMean::Finalize(AthenaArray<Real> const &w) {
  MeshBlock *pmb = pmy_block_;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  auto pthermo = Thermodynamics::GetInstance();

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
          for (int i = is; i <= ie; ++i) data(n, k, j, i) = w(n, k, j, i);

    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i)
          data(IPR, k, j, i) = pthermo->GetTemp(pmb, k, j, i);
  }

  // clear cycle
  ncycle_ = 0;
}
