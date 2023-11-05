//! \file particle_restart.cpp
//! \brief Implements functions derived from clas RestartGroup

// Athena++
#include <athena/mesh/mesh.hpp>

// canoe
#include <configure.hpp>

// n-body
#include "particle_data.hpp"
#include "particles.hpp"

// mpi
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

size_t ParticleBase::RestartDataSizeInBytes(Mesh const *pm) const {
  size_t size = sizeof(int);
  size += pc.size() * sizeof(ParticleData);

  // gather maximum size across all local blocks
  for (int b = 0; b < pm->nbtotal; b++) {
    auto pmb = pm->my_blocks(b);
    for (auto &pt : pmb->pimpl->all_particles) {
      size_t tmp = sizeof(int) + pt->pc.size() * sizeof(ParticleData);
      if (tmp > size) size = tmp;
    }
  }

  // gather maximum size across all MPI ranks
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, &size, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif

  return size;
}

//! \todo mark out invalid particles
void ParticleBase::DumpRestartData(char *pdst) const {
  int size = pc.size();

  std::memcpy(pdst, &size, sizeof(int));
  pdst += sizeof(int);
  std::memcpy(pdst, pc.data(), size * sizeof(ParticleData));
  pdst += size * sizeof(ParticleData);
}

size_t ParticleBase::LoadRestartData(char *psrc) {
  int size;

  std::memcpy(&size, psrc, sizeof(int));
  psrc += sizeof(int);
  pc.resize(size);
  std::memcpy(pc.data(), psrc, size * sizeof(ParticleData));
  psrc += size * sizeof(ParticleData);

  return size;
}
