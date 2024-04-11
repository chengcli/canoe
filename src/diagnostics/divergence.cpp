// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/reconstruct/interpolation.hpp>

// canoe
#include "diagnostics.hpp"

Divergence::Divergence(MeshBlock *pmb) : Diagnostics(pmb, "div,div_h") {
  type = "VECTORS";

  data.NewAthenaArray(2, ncells3_, ncells2_, ncells1_);
  v1f1_.NewAthenaArray(ncells3_, ncells2_, ncells1_ + 1);

  if (pmb->block_size.nx2 > 1) {
    v2f2_.NewAthenaArray(ncells3_, ncells2_ + 1, ncells1_);
  }

  if (pmb->block_size.nx3 > 1) {
    v3f3_.NewAthenaArray(ncells3_ + 1, ncells2_, ncells1_);
  }
}

void Divergence::Finalize(MeshBlock *pmb) {
  Coordinates *pcoord = pmb->pcoord;
  auto const &w = pmb->phydro->w;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  // interpolate velocities to cell faces
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie + 1; ++i) {
        v1f1_(k, j, i) = interp_cp4(w(IVX, k, j, i - 2), w(IVX, k, j, i - 1),
                                    w(IVX, k, j, i), w(IVX, k, j, i + 1));
      }
  if (pmb->block_size.nx2 > 1)
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je + 1; ++j)
        for (int i = is; i <= ie; ++i) {
          v2f2_(k, j, i) = interp_cp4(w(IVY, k, j - 2, i), w(IVY, k, j - 1, i),
                                      w(IVY, k, j, i), w(IVY, k, j + 1, i));
        }
  if (pmb->block_size.nx3 > 1)
    for (int k = ks; k <= ke + 1; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          v3f3_(k, j, i) = interp_cp4(w(IVZ, k - 2, j, i), w(IVZ, k - 1, j, i),
                                      w(IVZ, k, j, i), w(IVZ, k + 1, j, i));
        }

  // divergence
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pcoord->CellVolume(k, j, is, ie, vol_);
      for (int i = is; i <= ie; ++i) {
        pcoord->Face1Area(k, j, is, ie + 1, x1area_);
        data(0, k, j, i) =
            x1area_(i + 1) * v1f1_(k, j, i + 1) - x1area_(i) * v1f1_(k, j, i);

        if (pmb->block_size.nx2 > 1) {
          pcoord->Face2Area(k, j, is, ie, x2area_);
          pcoord->Face2Area(k, j + 1, is, ie, x2area_p1_);
          data(0, k, j, i) +=
              x2area_p1_(i) * v2f2_(k, j + 1, i) - x2area_(i) * v2f2_(k, j, i);
          data(1, k, j, i) =
              x2area_p1_(i) * v2f2_(k, j + 1, i) - x2area_(i) * v2f2_(k, j, i);
        }

        if (pmb->block_size.nx3 > 1) {
          pcoord->Face3Area(k, j, is, ie, x3area_);
          pcoord->Face3Area(k + 1, j, is, ie, x3area_p1_);
          data(0, k, j, i) +=
              x3area_p1_(i) * v3f3_(k + 1, j, i) - x3area_(i) * v3f3_(k, j, i);
          data(1, k, j, i) +=
              x3area_p1_(i) * v3f3_(k + 1, j, i) - x3area_(i) * v3f3_(k, j, i);
        }

        data(0, k, j, i) /= vol_(i);
        data(1, k, j, i) /= vol_(i);
      }
    }
}
