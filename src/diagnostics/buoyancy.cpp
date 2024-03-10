// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/hydro/srcterms/hydro_srcterms.hpp>
#include <athena/reconstruct/interpolation.hpp>

// diagnostics
#include "diagnostics.hpp"

Buoyancy::Buoyancy(MeshBlock* pmb) : Diagnostics(pmb, "b") {
  type = "SCALARS";
  data.NewAthenaArray(ncells3_, ncells2_, ncells1_);
  pf_.NewAthenaArray(ncells3_, ncells2_, ncells1_ + 1);
  grav_ = pmb->phydro->hsrc.GetG1();
}

void Buoyancy::Finalize(MeshBlock* pmb) {
  Coordinates* pcoord = pmb->pcoord;
  auto const& w = pmb->phydro->w;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie + 1; ++i) {
        // pf_(k,j,i) =
        // interp_cp4(w(IPR,k,j,i-2),w(IPR,k,j,i-1),w(IPR,k,j,i),w(IPR,k,j,i+1));
        pf_(k, j, i) = sqrt(w(IPR, k, j, i - 1) * w(IPR, k, j, i));
      }

  // buoyancy
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real dx = pmb->pcoord->dx1f(i);
        data(k, j, i) =
            (pf_(k, j, i) - pf_(k, j, i + 1)) / (dx * w(IDN, k, j, i)) + grav_;
      }

  // fix boundary condition
  bool has_top_neighbor = false;
  bool has_bot_neighbor = false;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock& nb = pmb->pbval->neighbor[n];
    if ((nb.ni.ox1 == -1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0))
      has_bot_neighbor = true;
    if ((nb.ni.ox1 == 1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0))
      has_top_neighbor = true;
  }

  if (!has_bot_neighbor) {
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) data(k, j, is) = data(k, j, is + 1);
  }

  if (!has_top_neighbor) {
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) data(k, j, ie) = data(k, j, ie - 1);
  }
}
