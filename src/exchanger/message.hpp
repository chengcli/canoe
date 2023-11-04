#ifndef SRC_COMMUNICATOR_COMMUNICATOR_HPP_
#define SRC_COMMUNICATOR_COMMUNICATOR_HPP_

// C/C++
#include <string>
#include <vector>

// canoe
#include <common.hpp>
#include <configure.hpp>

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

class MeshBlock;
struct NeighborBlock;

//! \brief Trait class providing Message information for class T
template <typename T>
struct Message {
  using DataType = double;

  constexpr static int num_buffers = 56;
  constexpr static std::string name;

#ifdef MPI_PARALLEL
  static MPI_Datatype mpi_type;
#endif
};

namespace MessageHelper {

int mpi_tag_ub;
int create_mpi_tag(int lid, int tid, std::string name);

}  // namespace MessageHelper

#endif  // SRC_COMMUNICATOR_COMMUNICATOR_HPP_
