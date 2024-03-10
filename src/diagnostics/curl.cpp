// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/reconstruct/interpolation.hpp>

// canoe
#include "diagnostics.hpp"

Curl::Curl(MeshBlock *pmb) : Diagnostics(pmb, "curl") {
  if (pmb->block_size.nx3 > 1) {
    type = "VECTORS";
    data.NewAthenaArray(3, ncells3_, ncells2_, ncells1_);
    v3f2_.NewAthenaArray(ncells3_ + 1, ncells2_ + 1, ncells1_ + 1);
    v2f3_.NewAthenaArray(ncells3_ + 1, ncells2_ + 1, ncells1_ + 1);
    v1f3_.NewAthenaArray(ncells3_ + 1, ncells2_ + 1, ncells1_ + 1);
    v3f1_.NewAthenaArray(ncells3_ + 1, ncells2_ + 1, ncells1_ + 1);
    v2f1_.NewAthenaArray(ncells3_ + 1, ncells2_ + 1, ncells1_ + 1);
    v1f2_.NewAthenaArray(ncells3_ + 1, ncells2_ + 1, ncells1_ + 1);
  } else if (pmb->block_size.nx2 > 1) {
    type = "SCALARS";
    data.NewAthenaArray(ncells3_, ncells2_, ncells1_);
    v2f1_.NewAthenaArray(ncells3_, ncells2_ + 1, ncells1_ + 1);
    v1f2_.NewAthenaArray(ncells3_, ncells2_ + 1, ncells1_ + 1);
  }
}

void Curl::Finalize(MeshBlock *pmb) {
  Coordinates *pcoord = pmb->pcoord;
  auto const &w = pmb->phydro->w;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  // interpolate velocities to cell faces
  if (pmb->block_size.nx3 > 1) {
    for (int k = ks; k <= ke + 1; ++k)
      for (int j = js; j <= je + 1; ++j)
        for (int i = is; i <= ie + 1; ++i) {
          v3f2_(k, j, i) = interp_cp4(w(IVZ, k, j - 2, i), w(IVZ, k, j - 1, i),
                                      w(IVZ, k, j, i), w(IVZ, k, j + 1, i));
          v2f3_(k, j, i) = interp_cp4(w(IVY, k - 2, j, i), w(IVY, k - 1, j, i),
                                      w(IVY, k, j, i), w(IVY, k + 1, j, i));
          v1f3_(k, j, i) = interp_cp4(w(IVX, k - 2, j, i), w(IVX, k - 1, j, i),
                                      w(IVX, k, j, i), w(IVX, k + 1, j, i));
          v3f1_(k, j, i) = interp_cp4(w(IVZ, k, j, i - 2), w(IVZ, k, j, i - 1),
                                      w(IVZ, k, j, i), w(IVZ, k, j, i + 1));
          v2f1_(k, j, i) = interp_cp4(w(IVY, k, j, i - 2), w(IVY, k, j, i - 1),
                                      w(IVY, k, j, i), w(IVY, k, j, i + 1));
          v1f2_(k, j, i) = interp_cp4(w(IVX, k, j - 2, i), w(IVX, k, j - 1, i),
                                      w(IVX, k, j, i), w(IVX, k, j + 1, i));
        }
  } else if (pmb->block_size.nx2 > 1) {
    for (int j = js; j <= je + 1; ++j)
      for (int i = is; i <= ie + 1; ++i) {
        v2f1_(ks, j, i) = interp_cp4(w(IVY, ks, j, i - 2), w(IVY, ks, j, i - 1),
                                     w(IVY, ks, j, i), w(IVY, ks, j, i + 1));
        v1f2_(ks, j, i) = interp_cp4(w(IVX, ks, j - 2, i), w(IVX, ks, j - 1, i),
                                     w(IVX, ks, j, i), w(IVX, ks, j + 1, i));
      }
  }

  // curl, FIXME, check whether this works for other coordinates
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      if (pmb->block_size.nx3 > 1) {
        pcoord->Edge2Length(k, j, is, ie, x2edge_);
        pcoord->Edge2Length(k + 1, j, is, ie, x2edge_p1_);
        pcoord->Edge3Length(k, j, is, ie, x3edge_);
        pcoord->Edge3Length(k, j + 1, is, ie, x3edge_p1_);
        pcoord->VolCenterFace1Area(k, j, is, ie, x1area_);

        for (int i = is; i <= ie; ++i) {
          data(0, k, j, i) = -(x3edge_p1_(i) * v3f2_(k, j + 1, i) -
                               x3edge_(i) * v3f2_(k, j, i) -
                               x2edge_p1_(i) * v2f3_(k + 1, j, i) +
                               x2edge_(i) * v2f3_(k, j, i)) /
                             x1area_(i);
        }
      }

      // curl
      if (pmb->block_size.nx3 > 1) {
        pcoord->Edge3Length(k, j, is, ie + 1, x3edge_);
        pcoord->Edge1Length(k, j, is, ie, x1edge_);
        pcoord->Edge1Length(k + 1, j, is, ie, x1edge_p1_);
        pcoord->VolCenterFace2Area(k, j, is, ie, x2area_);
        for (int i = is; i <= ie; ++i)
          data(1, k, j, i) = -(x1edge_p1_(i) * v1f3_(k + 1, j, i) -
                               x1edge_(i) * v1f3_(k, j, i) -
                               x3edge_(i + 1) * v3f1_(k, j, i + 1) +
                               x3edge_(i) * v3f1_(k, j, i)) /
                             x2area_(i);
      }

      // curl
      if (pmb->block_size.nx2 > 1) {
        pcoord->Edge1Length(k, j, is, ie, x1edge_);
        pcoord->Edge1Length(k, j + 1, is, ie, x1edge_p1_);
        pcoord->Edge2Length(k, j, is, ie + 1, x2edge_);
        pcoord->VolCenterFace3Area(k, j, is, ie, x3area_);
        if (pmb->block_size.nx3 > 1) {  // curl is a vector
          for (int i = is; i <= ie; ++i)
            data(2, k, j, i) = -(x2edge_(i + 1) * v2f1_(k, j, i + 1) -
                                 x2edge_(i) * v2f1_(k, j, i) -
                                 x1edge_p1_(i) * v1f2_(k, j + 1, i) +
                                 x1edge_(i) * v1f2_(k, j, i)) /
                               x3area_(i);
        } else {  // curl is a scalar
          for (int i = is; i <= ie; ++i)
            data(k, j, i) = -(x2edge_(i + 1) * v2f1_(k, j, i + 1) -
                              x2edge_(i) * v2f1_(k, j, i) -
                              x1edge_p1_(i) * v1f2_(k, j + 1, i) +
                              x1edge_(i) * v1f2_(k, j, i)) /
                            x3area_(i);
        }
      }
    }
}
