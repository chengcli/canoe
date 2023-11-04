#ifndef SRC_NBODY_PARTICLE_DATA_HPP_
#define SRC_NBODY_PARTICLE_DATA_HPP_

// C/C++ headers
#include <array>
#include <iosfwd>
#include <vector>

// Athena++ headers
#include <athena/athena.hpp>

// canoe
#include <configure.hpp>

// communication
#include <communication/neighbor_exchanger.hpp>

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif  // MPI_PARALLEL

class MeshBlock;

struct ParticleData {
  //! particle can form a linked list
  ParticleData* next;

  //! particle id, type id
  int pid, tid;

  //! extra integer data
  std::array<int, NINT_PARTICLE_DATA> ii;

  //! time instantiated, weight and charge
  Real time, weight, charge;

  //! positions
  Real x1, x2, x3;

  //! velocities
  Real v1, v2, v3;

  //! accelerations
  Real a1, a2, a3;

  //! extra real data
  std::array<Real, NREAL_PARTICLE_DATA> rr;
};

std::ostream& operator<<(std::ostream& os, ParticleData const& mp);

// Specialization for ParticleData Exchanger
template <>
struct NeighborExchanger<ParticleData> {
  using container_t = std::vector<ParticleData>;

  constexpr static int num_buffers = 56;
  constexpr static std::string name = "ParticleData";

#ifdef MPI_PARALLEL
  static MPI_Datatype mpi_type;
#endif  // MPI_PARALLEL
};

// helper functions
bool check_in_meshblock(ParticleData const& pd, MeshBlock const* pmb);
void commit_mpi_particle_data();
void free_mpi_particle_data();

#endif  // SRC_NBODY_PARTICLE_DATA_HPP_
