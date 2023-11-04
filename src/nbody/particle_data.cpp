//! \file particle_data.cpp
//! \brief Implementation of functions in class ParticleData

// C/C++ headers
#include <iostream>

// Athena++
#include <athena/mesh/mesh.hpp>

// canoe
#include <configure.hpp>

// exchanger
#include <exchanger/message_traits.hpp>

// nbody
#include "particle_data.hpp"
#include "particles.hpp"

std::ostream& operator<<(std::ostream& os, ParticleData const& pt) {
  os << "pid: " << pt.pid << ", tid: " << pt.tid << std::endl
     << "time: " << pt.time << ", weight: " << pt.weight
     << ", charge: " << pt.charge << std::endl
     << "x1: " << pt.x1 << ", v1: " << pt.v1 << ", a1: " << pt.a1 << std::endl
     << "x2: " << pt.x2 << ", v2: " << pt.v2 << ", a2: " << pt.a2 << std::endl
     << "x3: " << pt.x3 << ", v3: " << pt.v3 << ", a3: " << pt.a3 << std::endl;
  os << "extra real data: ";
  for (int i = 0; i < NREAL_PARTICLE_DATA; ++i) os << pt.rr[i] << ", ";
  os << std::endl;

  os << "extra integer data: ";
  for (int i = 0; i < NREAL_PARTICLE_DATA; ++i) os << pt.ii[i] << ", ";
  os << std::endl;

  return os;
}

namespace ParticlesHelper {

bool check_in_meshblock(ParticleData const& pd, MeshBlock const* pmb) {
  auto pm = pmb->pmy_mesh;

  bool inblock_x1 =
      (pd.x1 > pmb->block_size.x1min) && (pd.x1 < pmb->block_size.x1max);
  bool inblock_x2 =
      (pd.x2 > pmb->block_size.x2min) && (pd.x2 < pmb->block_size.x2max);
  bool inblock_x3 =
      (pd.x3 > pmb->block_size.x3min) && (pd.x3 < pmb->block_size.x3max);

  return inblock_x1 && (!pm->f2 || inblock_x2) && (!pm->f3 || inblock_x3);
}

#ifdef MPI_PARALLEL
#include <mpi.h>

void commit_mpi_particle_data() {
  int counts[3] = {1, 2 + NINT_PARTICLE_DATA, 8 + NREAL_PARTICLE_DATA};
  MPI_Datatype types[3] = {MPI_AINT, MPI_INT, MPI_ATHENA_REAL};
  MPI_Aint disps[3] = {offsetof(ParticleData, next),
                       offsetof(ParticleData, pid),
                       offsetof(ParticleData, time)};

  MPI_Type_create_struct(3, counts, disps, types,
                         &MessageTraits<ParticleBase>::mpi_type);
  MPI_Type_commit(&MessageTraits<ParticleBase>::mpi_type);
}

void free_mpi_particle_data() {
  MPI_Type_free(&MessageTraits<ParticleBase>::mpi_type);
}

#else  // NOT_MPI_PARALLEL

void commit_mpi_particle_data() {}
void free_mpi_particle_data() {}

#endif  // MPI_PARALLEL

}  // namespace ParticlesHelper
