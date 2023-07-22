// C/C++ headers
#include <sstream>
// #include <cstddef>
#include <functional>  // hash
#include <iostream>

// Athena++ classes headers
#include "../bvals/bvals.hpp"
#include "../debugger/debugger.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"
#include "material_point.hpp"
#include "particle_buffer.hpp"
#include "particles.hpp"

#ifdef MPI_PARALLEL
// defined in main.cpp
extern MPI_Datatype MPI_PARTICLE;
#endif

ParticleBuffer::ParticleBuffer(Particles *ppart) : pmy_particle(ppart) {
  for (int i = 0; i < 56; ++i) {
    particle_flag_[i] = BoundaryStatus::waiting;

#ifdef MPI_PARALLEL
    req_particle_send_[i] = MPI_REQUEST_NULL;
    req_particle_recv_[i] = MPI_REQUEST_NULL;
#endif
  }
}

ParticleBuffer::~ParticleBuffer() {}

int ParticleBuffer::CreateMPITag(int lid, int tid) {
  int TAG_PARTICLE = 15;
  int tag = BoundaryBase::CreateBvalsMPITag(lid, tid, TAG_PARTICLE);
  std::string str = pmy_particle->myname + std::to_string(tag);
  return std::hash<std::string>{}(str) % (Globals::mpi_tag_ub);
}

void ParticleBuffer::SendParticle() {
  Debugger *pdebug = pmy_particle->pmy_block->pdebug;
  for (int n = 0; n < pmy_particle->pmy_block->pbval->nneighbor; ++n) {
    NeighborBlock &nb = pmy_particle->pmy_block->pbval->neighbor[n];

    if (nb.snb.rank == Globals::my_rank) {  // on the same process
      MeshBlock *pmb =
          pmy_particle->pmy_block->pmy_mesh->FindMeshBlock(nb.snb.gid);
      pmb->ppart->ppb->particle_recv_[nb.targetid] = particle_send_[nb.bufid];
      pmb->ppart->ppb->particle_flag_[nb.targetid] = BoundaryStatus::arrived;
    }
#ifdef MPI_PARALLEL
    else {  // MPI
      int tag = CreateMPITag(nb.snb.lid, nb.targetid);
      int ssize = particle_send_[nb.bufid].size();
      MPI_Isend(particle_send_[nb.bufid].data(), ssize, MPI_PARTICLE,
                nb.snb.rank, tag, MPI_COMM_WORLD,
                &req_particle_send_[nb.bufid]);
#if DEBUG_LEVEL > 3
      if (ssize > 0) {
        std::cout << "- block " << Globals::my_rank << " send " << ssize << " "
                  << pmy_particle->myname << std::endl;
      }
#endif
    }
#endif
  }
}

void ParticleBuffer::RecvParticle() {
  int rsize, tag;
  Debugger *pdebug = pmy_particle->pmy_block->pdebug;

  MeshBlock *pmb = pmy_particle->pmy_block;
#ifdef MPI_PARALLEL
  MPI_Status status;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock &nb = pmb->pbval->neighbor[n];
    if (nb.snb.rank == Globals::my_rank) continue;  // local boundary received

    if (particle_flag_[nb.bufid] == BoundaryStatus::waiting) {
      int tag = CreateMPITag(pmb->lid, nb.bufid);
      MPI_Probe(nb.snb.rank, tag, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, MPI_PARTICLE, &rsize);
      particle_recv_[nb.bufid].resize(rsize);
      MPI_Irecv(particle_recv_[nb.bufid].data(), rsize, MPI_PARTICLE,
                nb.snb.rank, tag, MPI_COMM_WORLD,
                &req_particle_recv_[nb.bufid]);
#if DEBUG_LEVEL > 3
      if (rsize > 0) {
        std::cout << "- block " << Globals::my_rank << " receive " << rsize
                  << " " << pmy_particle->myname << std::endl;
      }
#endif
    }
  }
#endif
}

void ParticleBuffer::ClearBoundary() {
  MeshBlock *pmb = pmy_particle->pmy_block;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock &nb = pmb->pbval->neighbor[n];
    particle_flag_[nb.bufid] = BoundaryStatus::waiting;
    particle_send_[nb.bufid].clear();
    particle_recv_[nb.bufid].clear();

#ifdef MPI_PARALLEL
    if (nb.snb.rank != Globals::my_rank)
      MPI_Wait(&req_particle_send_[nb.bufid], MPI_STATUS_IGNORE);
#endif
  }
}
