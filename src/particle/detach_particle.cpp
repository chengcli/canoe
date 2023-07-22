/** @file detach_particle.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Monday May 31, 2021 13:38:00 PDT
 * @bug No known bugs.
 */

// C/C++ headers
// #include <algorithm>
// #include <iostream>

// Athena++ headers
#include "../bvals/bvals.hpp"
#include "../debugger/debugger.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"
#include "particle_buffer.hpp"
#include "particles.hpp"

// This subroutine will remove inactive particles (id < 0) and move particles to
// appropriate buffers if they moved out of the current domain
void ParticleBuffer::DetachParticle(std::vector<MaterialPoint> &mp) {
  MeshBlock *pmb = pmy_particle->pmy_block;
#if DEBUG_LEVEL > 2
  pmb->pdebug->Call("ParticleBuffer::DetachParticle-" + pmy_particle->myname);
  pmb->pdebug->CheckParticleConservation(pmy_particle->cnames_, mp);
#endif
  Mesh *pm = pmb->pmy_mesh;

  Real x1min = pmb->block_size.x1min;
  Real x1max = pmb->block_size.x1max;
  Real x2min = pmb->block_size.x2min;
  Real x2max = pmb->block_size.x2max;
  Real x3min = pmb->block_size.x3min;
  Real x3max = pmb->block_size.x3max;

  int ox1 = 0, ox2 = 0, ox3 = 0, fi1 = 0, fi2 = 0;
  std::vector<MaterialPoint>::iterator qi = mp.begin();
  std::vector<MaterialPoint>::iterator qj = mp.end();

  while (qi < qj) {
    // if particle is inactive, swap the current one with the last one
    if (qi->id < 0) {
      *qi = *(qj - 1);
      qj--;
      continue;
    }  // proceed to living particle
    // take care of reflective boundary condition
    if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::reflect &&
        qi->x1 < x1min) {
      qi->x1 = 2 * x1min - qi->x1;
      qi->v1 = -qi->v1;
    }
    if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::reflect &&
        qi->x1 > x1max) {
      qi->x1 = 2 * x1max - qi->x1;
      qi->v1 = -qi->v1;
    }
    ox1 = qi->x1 < x1min ? -1 : (qi->x1 > x1max ? 1 : 0);

    if (pm->f2) {
      if (pmb->pbval->block_bcs[inner_x2] == BoundaryFlag::reflect &&
          qi->x2 < x2min) {
        qi->x2 = 2 * x2min - qi->x2;
        qi->v2 = -qi->v2;
      } else if (pmb->pbval->block_bcs[inner_x2] == BoundaryFlag::polar &&
                 qi->x2 < 0.) {
        // \todo TODO: fix pole problem
      }
      if (pmb->pbval->block_bcs[outer_x2] == BoundaryFlag::reflect &&
          qi->x2 > x2max) {
        qi->x2 = 2 * x2max - qi->x2;
        qi->v2 = -qi->v2;
      } else if (pmb->pbval->block_bcs[outer_x2] == BoundaryFlag::polar &&
                 qi->x2 > 2. * M_PI) {
        // \todo TODO: fix pole problem
      }
      ox2 = qi->x2 < x2min ? -1 : (qi->x2 > x2max ? 1 : 0);
    }

    if (pm->f3) {
      if (pmb->pbval->block_bcs[inner_x3] == BoundaryFlag::reflect &&
          qi->x3 < x3min) {
        qi->x3 = 2 * x3min - qi->x3;
        qi->v3 = -qi->v3;
      }
      if (pmb->pbval->block_bcs[outer_x3] == BoundaryFlag::reflect &&
          qi->x3 > x3max) {
        qi->x3 = 2 * x3max - qi->x3;
        qi->v3 = -qi->v3;
      }
      ox3 = qi->x3 < x3min ? -1 : (qi->x3 > x3max ? 1 : 0);
    }

    if (pm->multilevel) {
      // reserved implementation for multilevel, fi1, fi2
    }

    int bid = BoundaryBase::FindBufferID(ox1, ox2, ox3, fi1, fi2);
    if (bid == -1) {  // particle inside domain
      qi++;
    } else {  // particle moved out of the domain
      // std::cout << "particle (" << qi->id << "," << qi->x2 << ") move to " <<
      // bid << std::endl;; std::cout << "block range = (" << x1min << "," <<
      // x1max << ")" << std::endl;
      particle_send_[bid].push_back(*qi);
      // pmb->pdebug->msg << *qi;
      *qi = *(qj - 1);
      qj--;
    }
  }

  // particles beyond qi are inactive particles. Remove them from the list
  mp.resize(qi - mp.begin());
#if DEBUG_LEVEL > 2
  pmb->pdebug->CheckParticleConservation(pmy_particle->cnames_, mp);
  pmb->pdebug->Leave();
#endif
}
