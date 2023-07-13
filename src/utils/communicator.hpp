/** \file communicator.hpp
 * \brief
 *
 * \author Cheng Li (chengcli@umich.edu)
 * \date Thursday Apr 14, 2022 11:39:35 EDT
 */

#ifndef SRC_UTILS_COMMUNICATOR_HPP
#define SRC_UTILS_COMMUNICATOR_HPP

// athena
#include <athena/athena.hpp>

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

class MeshBlock;
struct NeighborBlock;

class Communicator {
 public:
  Communicator(MeshBlock *pmb);
  ~Communicator();
  int getRank(CoordinateDirection dir) const;
  void setColor(CoordinateDirection dir);
  // void reduceData23(Real *send, Real *recv);
  void gatherData(Real *send, Real *recv, int size) const;
  void gatherDataInPlace(Real *recv, int size) const;

  NeighborBlock const *findBotNeighbor() const;
  NeighborBlock const *findTopNeighbor() const;
  NeighborBlock const *findLeftNeighbor() const;
  NeighborBlock const *findRightNeighbor() const;
  NeighborBlock const *findBackNeighbor() const;
  NeighborBlock const *findFrontNeighbor() const;

 private:
  MeshBlock *pmy_block_;
  int *color_;
  int *brank_;
#ifdef MPI_PARALLEL
  MPI_Comm comm_;
#endif
};

#endif  // SRC_UTILS_COMMUNICATOR_HPP
