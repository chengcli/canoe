#ifndef SRC_EXCHANGER_MESSAGE_TRAITS_HPP_
#define SRC_EXCHANGER_MESSAGE_TRAITS_HPP_

// C/C++
#include <string>

// canoe
#include <configure.hpp>

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

class RadiationBand;
class ParticleData;
class ParticleBase;

//! \brief Traits class providing Message information for class T
template <typename T>
struct MessageTraits {
  using DataType = double;

  constexpr static int num_buffers = 1;
  constexpr static const char* name = "";

#ifdef MPI_PARALLEL
  static MPI_Datatype mpi_type;
#endif
};

// Specialization for RadiationBand Exchanger
template <>
struct MessageTraits<RadiationBand> {
  using DataType = Real;

  constexpr static int num_buffers = 2;
  constexpr static const char* name = "RadiationBand";

#ifdef MPI_PARALLEL
  static MPI_Datatype mpi_type;
#endif  // MPI_PARALLEL
};

// Specialization for Particle Exchanger
template <>
struct MessageTraits<ParticleBase> {
  using DataType = ParticleData;

  constexpr static int num_buffers = 56;
  constexpr static const char* name = "ParticleBase";

#ifdef MPI_PARALLEL
  static MPI_Datatype mpi_type;
#endif  // MPI_PARALLEL
};

namespace MessageHelper {

extern int mpi_tag_ub;
int create_mpi_tag(int lid, int tid, std::string name);

}  // namespace MessageHelper

#endif  // SRC_COMMUNICATOR_COMMUNICATOR_HPP_
