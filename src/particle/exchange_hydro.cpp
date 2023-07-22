/** @file exchange_hydro.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Thursday Jun 17, 2021 07:04:02 PDT
 * @bug No known bugs.
 */

// Athena++ headers
#include "../debugger/debugger.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../math/interpolation.h"  // locate
#include "../mesh/mesh.hpp"
#include "particles.hpp"

void Particles::ExchangeHydro(std::vector<MaterialPoint> &mp,
                              AthenaArray<Real> &du, AthenaArray<Real> const &w,
                              Real dt) {
  MeshBlock *pmb = pmy_block;
  pmb->pdebug->Call("Particles::ExchangeHydro-" + myname);

  Mesh *pm = pmb->pmy_mesh;
  AthenaArray<Real> v1, v2, v3;
  Real loc[3];

  v1.InitWithShallowSlice(const_cast<AthenaArray<Real> &>(w), 4, IM1, 1);
  v2.InitWithShallowSlice(const_cast<AthenaArray<Real> &>(w), 4, IM2, 1);
  v3.InitWithShallowSlice(const_cast<AthenaArray<Real> &>(w), 4, IM3, 1);

  Real g1 = pmb->phydro->hsrc.GetG1();
  Real g2 = pmb->phydro->hsrc.GetG2();
  Real g3 = pmb->phydro->hsrc.GetG3();

  // loop over particles
  for (std::vector<MaterialPoint>::iterator q = mp.begin(); q != mp.end();
       ++q) {
    loc[0] = q->x3;
    loc[1] = q->x2;
    loc[2] = q->x1;

    int f = pm->f2 + pm->f3;
    // interpolate Eulerian velocity to particle velocity
    interpn(&q->v1, loc + 2 - f, v1.data(), xcenter_.data() + 2 - f,
            dims_.data() + 2 - f, 1 + f);
    if (pm->f2)
      interpn(&q->v2, loc + 2 - f, v2.data(), xcenter_.data() + 2 - f,
              dims_.data() + 2 - f, 1 + f);
    else
      q->v2 = 0.;
    if (pm->f3)
      interpn(&q->v3, loc + 2 - f, v3.data(), xcenter_.data() + 2 - f,
              dims_.data() + 2 - f, 1 + f);
    else
      q->v3 = 0.;

    assert(!std::isnan(q->v1));
    assert(!std::isnan(q->v2));
    assert(!std::isnan(q->v3));
  }
  pmb->pdebug->Leave();
}
