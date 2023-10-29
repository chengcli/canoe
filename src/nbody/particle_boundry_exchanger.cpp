//! \file particle_boundry_exchanger.cpp
//! \brief Implements functions derived from class BoundaryExchanger

// C/C++ header
#include <iostream>
#include <sstream>
#include <stdexcept>

// application
#include <application/exceptions.hpp>

// Athena++ headers
#include <athena/bvals/bvals.hpp>
#include <athena/mesh/mesh.hpp>

// n-body
#include "particle_data.hpp"
#include "particles.hpp"

bool ParticleBase::AttachTo(ParticleContainer &container) {
  bool success = true;
  int test;

  auto pmb = getMeshBlock();
  auto pm = pmb->pmy_mesh;

  for (auto &nb : pmb->pbval->neighbor) {
    if (status_flag_[nb.bufid] == BoundaryStatus::completed) continue;

#ifdef MPI_PARALLEL
    if (status_flag_[nb.bufid] == BoundaryStatus::waiting) {
      MPI_Test(&req_mpi_recv_[nb.bufid], &test, MPI_STATUS_IGNORE);
      if (test)
        status_flag_[nb.bufid] = BoundaryStatus::arrived;
      else {
        success = false;
        continue;
      }
    }
#endif

    if (status_flag_[nb.bufid] == BoundaryStatus::arrived) {
      auto it = recv_buffer_[nb.bufid].begin();

      for (; it != recv_buffer_[nb.bufid].end(); ++it) {
        // 0:INNER_X1, 1:OUTER_X1
        if (pm->mesh_bcs[(nb.ni.ox1 + 1) >> 1] == BoundaryFlag::periodic) {
          if (it->x1 >= pm->mesh_size.x1max)
            it->x1 -= pm->mesh_size.x1max - pm->mesh_size.x1min;
          if (it->x1 <= pm->mesh_size.x1min)
            it->x1 += pm->mesh_size.x1max - pm->mesh_size.x1min;
        }

        // 2:INNER_X2, 3:OUTER_X2
        if (pm->mesh_bcs[2 + ((nb.ni.ox2 + 1) >> 1)] == BoundaryFlag::periodic) {
          if (it->x2 >= pm->mesh_size.x2max)
            it->x2 -= pm->mesh_size.x2max - pm->mesh_size.x2min;
          if (it->x2 <= pm->mesh_size.x2min)
            it->x2 += pm->mesh_size.x2max - pm->mesh_size.x2min;
        }

        // 4:INNER_X3, 5:OUTER_X3
        if (pm->mesh_bcs[4 + ((nb.ni.ox3 + 1) >> 1)] == BoundaryFlag::periodic) {
          if (it->x3 >= pm->mesh_size.x3max)
            it->x3 -= pm->mesh_size.x3max - pm->mesh_size.x3min;
          if (it->x3 <= pm->mesh_size.x3min)
            it->x3 += pm->mesh_size.x3max - pm->mesh_size.x3min;
        }

        bool in_meshblock = check_in_meshblock(*it, pmb);
        if (!in_meshblock) {
          throw RuntimeError("ParticleBase::AttachParticles",
                             "Particle moved out of MeshBlock limits");
        } else {
          container.push_back(*it);
        }
      }

      SetBoundaryStatus(nb.bufid, BoundaryStatus::completed);  // completed
    }
  }

  return success;
}

// This subroutine will remove inactive particles (id < 0) and move particles to
// appropriate buffers if they moved out of the current domain
void ParticleBase::DetachTo(ParticleContainer &buffer) {
  auto pmb = getMeshBlock();
  auto pm = pmb->pmy_mesh;

  Real x1min = pmb->block_size.x1min;
  Real x1max = pmb->block_size.x1max;
  Real x2min = pmb->block_size.x2min;
  Real x2max = pmb->block_size.x2max;
  Real x3min = pmb->block_size.x3min;
  Real x3max = pmb->block_size.x3max;

  int ox1 = 0, ox2 = 0, ox3 = 0, fi1 = 0, fi2 = 0;
  auto qi = buffer.begin();
  auto qj = buffer.end();

  while (qi < qj) {
    // if particle is inactive, swap the current one with the last one
    if (qi->pid < 0) {
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
      send_buffer_[bid].push_back(*qi);
      // pmb->pdebug->msg << *qi;
      *qi = *(qj - 1);
      qj--;
    }
  }

  // particles beyond qi are inactive particles. Remove them from the list
  buffer.resize(qi - buffer.begin());
}
