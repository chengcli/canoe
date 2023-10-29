#ifndef SRC_COMMUNICATOR_COMMUNICATOR_HPP_
#define SRC_COMMUNICATOR_COMMUNICATOR_HPP_

// C/C++
#include <vector>

// canoe
#include <configure.hpp>

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

template <typename T>
struct Communicator {
  using container_t = std::vector<T>;

  constexpr static int uid = -1;
  constexpr static int max_neighbors = 56;

#ifdef MPI_PARALLEL
  static MPI_Datatype mpi_type;
  static MPI_Comm mpi_comm;
#endif
};

#ifdef MPI_PARALLEL
template <typename T>
MPI_Comm Communicator<T>::mpi_comm = MPI_COMM_WORLD;
#endif

#endif  // SRC_COMMUNICATOR_COMMUNICATOR_HPP_
