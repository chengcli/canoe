/** @file attach_particle.cpp
 * @brief attach particle to the particle chain
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Friday Jun 04, 2021 20:26:36 PDT
 * @bug No known bugs.
 */

// C/C++ header
#include <iostream>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../bvals/bvals.hpp"
#include "../debugger/debugger.hpp"
#include "../mesh/mesh.hpp"
#include "particle_buffer.hpp"
#include "particles.hpp"

bool ParticleBuffer::AttachParticle(std::vector<MaterialPoint> &mp) {
  bool success = true;
  int test;

  MeshBlock *pmb = pmy_particle->pmy_block;
#if DEBUG_LEVEL > 2
  pmb->pdebug->Call("ParticleBuffer::AttachParticle-" + pmy_particle->myname);
#endif
  std::stringstream &msg = pmb->pdebug->msg;
  Mesh *pm = pmb->pmy_mesh;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock &nb = pmb->pbval->neighbor[n];
    if (particle_flag_[nb.bufid] == BoundaryStatus::completed) continue;

#ifdef MPI_PARALLEL
    if (particle_flag_[nb.bufid] == BoundaryStatus::waiting) {
      MPI_Test(&req_particle_recv_[nb.bufid], &test, MPI_STATUS_IGNORE);
      if (test)
        particle_flag_[nb.bufid] = BoundaryStatus::arrived;
      else {
        success = false;
        continue;
      }
    }
#endif

    if (particle_flag_[nb.bufid] == BoundaryStatus::arrived) {
      std::vector<MaterialPoint>::iterator it =
          particle_recv_[nb.bufid].begin();
#if DEBUG_LEVEL > 3
      if (particle_recv_[nb.bufid].size() > 0)
        std::cout << "- block " << Globals::my_rank << " attach "
                  << particle_recv_[nb.bufid].size() << std::endl;
#endif
      for (; it != particle_recv_[nb.bufid].end(); ++it) {
        // 0:INNER_X1, 1:OUTER_X1
        if (pm->mesh_bcs[nb.ni.ox1 + 1 >> 1] == BoundaryFlag::periodic) {
          if (it->x1 >= pm->mesh_size.x1max)
            it->x1 -= pm->mesh_size.x1max - pm->mesh_size.x1min;
          if (it->x1 <= pm->mesh_size.x1min)
            it->x1 += pm->mesh_size.x1max - pm->mesh_size.x1min;
        }

        // 2:INNER_X2, 3:OUTER_X2
        if (pm->mesh_bcs[2 + (nb.ni.ox2 + 1 >> 1)] == BoundaryFlag::periodic) {
          if (it->x2 >= pm->mesh_size.x2max)
            it->x2 -= pm->mesh_size.x2max - pm->mesh_size.x2min;
          if (it->x2 <= pm->mesh_size.x2min)
            it->x2 += pm->mesh_size.x2max - pm->mesh_size.x2min;
        }

        // 4:INNER_X3, 5:OUTER_X3
        if (pm->mesh_bcs[4 + (nb.ni.ox3 + 1 >> 1)] == BoundaryFlag::periodic) {
          if (it->x3 >= pm->mesh_size.x3max)
            it->x3 -= pm->mesh_size.x3max - pm->mesh_size.x3min;
          if (it->x3 <= pm->mesh_size.x3min)
            it->x3 += pm->mesh_size.x3max - pm->mesh_size.x3min;
        }

        bool out_of_domain = (it->x1 < pmb->block_size.x1min ||
                              it->x1 > pmb->block_size.x1max) ||
                             (pm->f2 && (it->x2 < pmb->block_size.x2min ||
                                         it->x2 > pmb->block_size.x2max)) ||
                             (pm->f3 && (it->x3 < pmb->block_size.x3min ||
                                         it->x3 > pmb->block_size.x3max));
        if (out_of_domain) {
// \todo once in while, a wire particle appears
// it has an id but density and positions are uninitialized
// I don't where it comes from and it will take time to debug
// leave it for now (solved, linked list problem in detach particle)
#if DEBUG_LEVEL > 3
          success = false;
          msg << "### FATAL ERROR in ParticleBuffer::AttachParticles. "
                 "Particles "
              << pmy_particle->myname << " moved out of MeshBlock limits"
              << std::endl;
          msg << (int)pm->mesh_bcs[0] << " " << (int)pm->mesh_bcs[1] << " "
              << (int)pm->mesh_bcs[2] << " " << (int)pm->mesh_bcs[3] << " "
              << (int)pm->mesh_bcs[4] << " " << (int)pm->mesh_bcs[5]
              << std::endl;
          msg << it->id << " " << it->type << " " << it->rho << std::endl;
          msg << it->x1 << " " << pmb->block_size.x1min << " "
              << pmb->block_size.x1max << std::endl;
          msg << it->x2 << " " << pmb->block_size.x2min << " "
              << pmb->block_size.x2max << std::endl;
          msg << it->x3 << " " << pmb->block_size.x3min << " "
              << pmb->block_size.x3max << std::endl;
          msg << it->v1 << std::endl;
          msg << it->v2 << std::endl;
          msg << it->v3 << std::endl;
          // for (auto x = particle_recv_[nb.bufid].begin(); x !=
          // particle_recv_[nb.bufid].end(); ++x)
          //   msg << x->id << " " << x->type << " " << x->rho << std::endl;
          ATHENA_ERROR(msg);
#endif
        } else {
          mp.push_back(*it);
          // msg << *it;
        }
      }
      particle_flag_[nb.bufid] = BoundaryStatus::completed;  // completed
    }
  }
#if DEBUG_LEVEL > 2
  pmb->pdebug->Leave();
#endif

  return success;
}
