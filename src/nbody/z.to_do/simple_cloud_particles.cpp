/** @file simple_cloud_particles.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Wednesday Jun 09, 2021 18:33:03 PDT
 * @bug No known bugs.
 */

// C/C++ header
#include <iostream>
#include <sstream>
#include <stdexcept>

// Athena++ header
#include "../coordinates/coordinates.hpp"
#include "../debugger/debugger.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../math/interpolation.h"  // locate
#include "../mesh/mesh.hpp"
#include "particles.hpp"

SimpleCloudParticles::SimpleCloudParticles(MeshBlock *pmb, ParameterInput *pin,
                                           std::string name)
    : Particles(pmb, pin, name, 2) {
  // ATHENA_LOG("SimpleCloudParticles");
  pmb->pdebug->Enter("SimpleCloudParticles");
  std::stringstream &msg = pmb->pdebug->msg;
  msg << "- first category is " << name + " cloud" << std::endl
      << "- second category is " << name + " precipitation" << std::endl;
  cnames_.resize(2);
  cnames_[0] = "cloud";
  cnames_[1] = "rain";

  Real mu = pin->GetReal("particles", name + ".mu");
  Real cc = pin->GetReal("particles", name + ".cc");

  mu_.push_back(mu);
  mu_.push_back(mu);

  cc_.push_back(cc);
  cc_.push_back(cc);
  pmb->pdebug->Leave();
}

void SimpleCloudParticles::ExchangeHydro(std::vector<MaterialPoint> &mp,
                                         AthenaArray<Real> &du,
                                         AthenaArray<Real> const &w, Real dt) {
  MeshBlock *pmb = pmy_block;
  pmb->pdebug->Call("SimpleCloudParticles::ExchangeHydro-" + myname);
  std::stringstream &msg = pmb->pdebug->msg;

  Mesh *pm = pmb->pmy_mesh;
  Coordinates *pcoord = pmb->pcoord;
  AthenaArray<Real> v1, v2, v3;
  Real loc[3];

  v1.InitWithShallowSlice(const_cast<AthenaArray<Real> &>(w), 4, IM1, 1);
  v2.InitWithShallowSlice(const_cast<AthenaArray<Real> &>(w), 4, IM2, 1);
  v3.InitWithShallowSlice(const_cast<AthenaArray<Real> &>(w), 4, IM3, 1);

  Real g1 = pmb->phydro->hsrc.GetG1();
  Real g2 = pmb->phydro->hsrc.GetG2();
  Real g3 = pmb->phydro->hsrc.GetG3();

  int f = pm->f2 + pm->f3;
  for (std::vector<MaterialPoint>::iterator q = mp.begin(); q != mp.end();
       ++q) {
    loc[0] = q->x3;
    loc[1] = q->x2;
    loc[2] = q->x1;

    interpnf(&q->v1, loc + 2 - f, v1.data(), xcenter_.data() + 2 - f,
             dims_.data() + 2 - f, 1 + f);
    if (pm->f2)
      interpnf(&q->v2, loc + 2 - f, v2.data(), xcenter_.data() + 2 - f,
               dims_.data() + 2 - f, 1 + f);
    else
      q->v2 = 0.;
    if (pm->f3)
      interpnf(&q->v3, loc + 2 - f, v3.data(), xcenter_.data() + 2 - f,
               dims_.data() + 2 - f, 1 + f);
    else
      q->v3 = 0.;

    if (std::isnan(q->v1) || std::isnan(q->v2) || std::isnan(q->v3)) {
      msg << "### FATAL ERROR in SimpleCloudParticle::ExchangeHydro. Particles "
          << myname << " velocity is nan" << std::endl;
      msg << "id = " << q->id << std::endl
          << "type = " << q->type << std::endl
          << "x1 = " << q->x1 << std::endl
          << "x2 = " << q->x2 << std::endl
          << "x3 = " << q->x3 << std::endl
          << "v1 = " << q->v1 << std::endl
          << "v2 = " << q->v2 << std::endl
          << "v3 = " << q->v3 << std::endl;
      Real k = locate(xcenter_.data(), q->x3, dims_[0]);
      Real j = locate(xcenter_.data() + dims_[0], q->x2, dims_[1]);
      Real i = locate(xcenter_.data() + dims_[0] + dims_[1], q->x1, dims_[2]);
      std::cout << i << " " << pcoord->x1v(i) << " " << pcoord->x1v(i + 1)
                << std::endl;
      std::cout << j << " " << pcoord->x2v(j) << " " << pcoord->x2v(j + 1)
                << std::endl;
      std::cout << k << " " << pcoord->x3v(k) << " " << pcoord->x3v(k + 1);
      ATHENA_ERROR(msg);
    }

    // \todo TODO: hard coded sedimentation
    if (q->type == 1) q->v1 += -10.;

    // add gravititional acceleration
    int k, j, i;
    k = locate(xface_.data(), q->x3, dims_[0] + 1);
    j = locate(xface_.data() + dims_[0] + 1, q->x2, dims_[1] + 1);
    i = locate(xface_.data() + dims_[0] + dims_[1] + 2, q->x1, dims_[2] + 1);

    Real src = dt * q->rho;
    du(IM1, k, j, i) += src * g1;
    du(IM2, k, j, i) += src * g2;
    du(IM3, k, j, i) += src * g3;
    du(IEN, k, j, i) += src * (g1 * q->v1 + g2 * q->v2 + g3 * q->v3);
  }

#if DEBUG_LEVEL > 2
  pmb->pdebug->CheckConservation("du", du, pmb->is, pmb->ie, pmb->js, pmb->je,
                                 pmb->ks, pmb->ke);
#endif
  pmb->pdebug->Leave();
}
