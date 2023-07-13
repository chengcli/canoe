/** \file communicator.hpp
 * \brief
 *
 * \author Cheng Li (chengcli@umich.edu)
 * \date Thursday Apr 14, 2022 11:39:35 EDT
 */

#ifndef SRC_UTILS_COMMUNICATOR_HPP
#define SRC_UTILS_COMMUNICATOR_HPP

// C/C++
#include <vector>

// athena
#include <athena/athena.hpp>

// canoe
#include <configure.hpp>

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

class MeshBlock;
struct NeighborBlock;

void find_neighbors(MeshBlock *pmb, CoordinateDirection dir,
                    NeighborBlock *bblock, NeighborBlock *tblock);

class Communicator {
 protected:
  //! Protected ctor access thru static member function Instance
  Communicator();

 public:
  ~Communicator();

  static Communicator *GetInstance();

  int GetRank(MeshBlock *pmb, CoordinateDirection dir) const;
  void SetColor(MeshBlock *pmb, CoordinateDirection dir);

  // void reduceData23(Real *send, Real *recv);
  void GatherData(Real *send, Real *recv, int size) const;
  void GatherDataInPlace(Real *recv, int size) const;

  NeighborBlock const *FindBotNeighbor(MeshBlock *pmb) const;
  NeighborBlock const *FindTopNeighbor(MeshBlock *pmb) const;
  NeighborBlock const *FindLeftNeighbor(MeshBlock *pmb) const;
  NeighborBlock const *FindRightNeighbor(MeshBlock *pmb) const;
  NeighborBlock const *FindBackNeighbor(MeshBlock *pmb) const;
  NeighborBlock const *FindFrontNeighbor(MeshBlock *pmb) const;

 private:
  std::vector<int> color_;
  std::vector<int> brank_;

  //! Pointer to the single Application instance
  static Communicator *mycomm_;

#ifdef MPI_PARALLEL
  MPI_Comm comm_;
#endif
};

#endif  // SRC_UTILS_COMMUNICATOR_HPP
