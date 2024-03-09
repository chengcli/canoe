// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/reconstruct/interpolation.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// canoe
#include "diagnostics.hpp"

HydroMean::HydroMean(MeshBlock *pmb) : Diagnostics(pmb, "mean") {
  type = "VECTORS";
  varname = "rho_bar,";
  long_name = "mean density,";
  units = "kg/m^3,";

  for (int n = 1; n <= NVAPOR; ++n) {
    varname += "q" + std::to_string(n) + "_bar,";
    units += "kg/kg,";
    long_name += "mean vapor,";
  }
  for (int n = 0; n < 3; ++n) {
    units += "m/s,";
    varname += "vel" + std::to_string(n + 1) + "_bar,";
    long_name += "mean velocity,";
  }
  units += "K";
  varname += "T_bar";
  long_name += "mean temperature";

  SetLongName(long_name);

  data.NewAthenaArray(NHYDRO, ncells3_, ncells2_, ncells1_);
}

HydroMean::~HydroMean() { data.DeleteAthenaArray(); }

void HydroMean::Progress(AthenaArray<Real> const &w) {
  MeshBlock *pmb = pmy_block_;

  if (ncycle == 0) std::fill(data.data(), data.data() + data.GetSize(), 0.);

  for (int n = 0; n < IPR; ++n)
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j)
        for (int i = pmb->is; i <= pmb->ie; ++i)
          data(n, k, j, i) += w(n, k, j, i);

  // Temperature field is save in IPR
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i)
        data(IPR, k, j, i) += pmb->pthermo->GetTemp(w.at(k, j, i));

  ncycle++;
}

void HydroMean::Finalize(AthenaArray<Real> const &w) {
  MeshBlock *pmb = pmy_block_;

  // hydro mean
  if (ncycle > 0) {
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = pmb->ks; k <= pmb->ke; ++k)
        for (int j = pmb->js; j <= pmb->je; ++j)
          for (int i = pmb->is; i <= pmb->ie; ++i) data(n, k, j, i) /= ncycle;
  } else {
    for (int n = 0; n < IPR; ++n)
      for (int k = pmb->ks; k <= pmb->ke; ++k)
        for (int j = pmb->js; j <= pmb->je; ++j)
          for (int i = pmb->is; i <= pmb->ie; ++i)
            data(n, k, j, i) = w(n, k, j, i);

    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j)
        for (int i = pmb->is; i <= pmb->ie; ++i)
          data(IPR, k, j, i) = pmb->pthermo->GetTemp(w.at(k, j, i));
  }

  // clear cycle
  ncycle = 0;
}
