// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// climath
#include <climath/core.h>

// canoe
#include <impl.hpp>

#include "parameterization.hpp"

namespace Parameterization {

Real top_sponge_lyr_width;
Real top_sponge_lyr_nu;

void init_parameterization(MeshBlock *pmb, ParameterInput *pin) {
  top_sponge_lyr_width =
      pin->GetOrAddReal("parameterization", "top_sponge_lyr.width", 0.);
  top_sponge_lyr_nu =
      pin->GetOrAddReal("parameterization", "top_sponge_lyr.nu", 0.);

  par_top_sponge_lyr_kl78(pmb);
}

void par_top_sponge_lyr_kl78(MeshBlock *pmb) {
  Real const &width = top_sponge_lyr_width;
  auto pcoord = pmb->pcoord;
  auto &nu = pmb->phydro->hdif.nu;

  Real x1max = pmb->block_size.x1max;
  Real ztop = x1max - width;
  for (int k = pmb->ks; k <= pmb->ke; k++)
    for (int j = pmb->js; j <= pmb->je; j++)
      for (int i = pmb->is; i <= pmb->ie; i++) {
        Real eta = (pcoord->x1f(i + 1) - ztop) / width;
        if (eta < 0) break;

        Real scale = sqr(sin(M_PI / 2. * eta));
        nu(0, k, j, i) = top_sponge_lyr_nu * scale;
      }
}

}  // namespace Parameterization
