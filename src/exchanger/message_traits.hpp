#ifndef SRC_EXCHANGER_MESSAGE_TRAITS_HPP_
#define SRC_EXCHANGER_MESSAGE_TRAITS_HPP_

// C/C++
#include <string>

// canoe
#include <configure.hpp>

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//! \brief Traits class providing Message information for class T
template <typename T>
struct MessageTraits {
  using DataType = double;

  constexpr static int num_buffers;
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
