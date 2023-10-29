/** @file aggregate_density.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Sunday Jun 13, 2021 12:48:16 PDT
 * @bug No known bugs.
 */

// C/C++ headers
#include <iostream>

// Athena++ headers
#include "../coordinates/coordinates.hpp"
#include "../debugger/debugger.hpp"
#include "../globals.hpp"
#include "../math/interpolation.h"  // locate
#include "../mesh/mesh.hpp"
#include "particles.hpp"

// AggregateDensity executes after AttachParticle and before chemistry
// mp stores all active particles (inactive particles are removed in
// DetachParticle) u stores the aggregated particle density and will be
// populated based on mp
void Particles::AggregateDensity(AthenaArray<Real> &u,
                                 std::vector<MaterialPoint> &mp) {
  MeshBlock *pmb = pmy_block;
  Coordinates *pcoord = pmb->pcoord;
  pmb->pdebug->Call("Particles::AggregateDensity-" + myname);

#if DEBUG_LEVEL > 2
  pmb->pdebug->CheckParticleConservation(cnames_, mp);
#endif

  // initialize cell aggregated particle density u and particle linked list
  // pcell_
  u.ZeroClear();
  std::fill(pcell_.data(), pcell_.data() + pcell_.GetSize(), nullptr);
  int i, j, k;

  // loop over particles
  for (std::vector<MaterialPoint>::iterator q = mp.begin(); q != mp.end();
       ++q) {
    k = locate(xface_.data(), q->x3, dims_[0] + 1);
    j = locate(xface_.data() + dims_[0] + 1, q->x2, dims_[1] + 1);
    i = locate(xface_.data() + dims_[0] + dims_[1] + 2, q->x1, dims_[2] + 1);

    // it may happen that the last cell is [a, b] and x = b
#if DEBUG_LEVEL > 3
    if (k < pmb->ks || k > pmb->ke || j < pmb->js || j > pmb->je ||
        i < pmb->is || i > pmb->ie) {
      std::stringstream msg;
      msg << "loc = (" << q->x1 << "," << q->x2 << "," << q->x3 << ")"
          << std::endl;
      << pcoord->x1f(i) << "," << pcoord->x1f(i + 1) << std::endl;
      << pcoord->x2f(j) << "," << pcoord->x2f(j + 1) << std::endl;
      << pcoord->x3f(k) << "," << pcoord->x3f(k + 1) << std::endl;
      ATHENA_ERROR(msg);
    }
#endif

    u(q->type, k, j, i) += q->rho;

    if (pcell_(q->type, k, j, i) == nullptr) {
      pcell_(q->type, k, j, i) = &(*q);
      pcell_(q->type, k, j, i)->next = nullptr;
    } else {  // insert sort from low to high density
      MaterialPoint *tmp;
      if (pcell_(q->type, k, j, i)->rho < q->rho) {
        MaterialPoint *pc = pcell_(q->type, k, j, i);
        while (pc->next != nullptr && pc->next->rho < q->rho) pc = pc->next;
        tmp = pc->next;
        pc->next = &(*q);
        pc->next->next = tmp;
      } else {
        tmp = pcell_(q->type, k, j, i);
        pcell_(q->type, k, j, i) = &(*q);
        pcell_(q->type, k, j, i)->next = tmp;
      }
    }
  }

#if DEBUG_LEVEL > 2
  pmb->pdebug->CheckConservation("u", u, pmb->is, pmb->ie, pmb->js, pmb->je,
                                 pmb->ks, pmb->ke);
#endif

  // copy u to u1_ such that u1_ represents the aggregated particle density
  // before chemistry for the next step
  u1_ = u;
  pmb->pdebug->Leave();
}
