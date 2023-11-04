#ifndef SRC_EXCHANGER_NEIGHBOR_EXCHANGER_HPP_
#define SRC_EXCHANGER_NEIGHBOR_EXCHANGER_HPP_

// athena
#include <athena/globals.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <configure.hpp>

// exchanger
#include "exchanger.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

template <typename T>
void NeighborExchanger<T>::Transfer(MeshBlock const *pmb, int n) {
  for (auto &nb : pmb->pbval->neighbor) {
    if (nb.snb.rank == Globals::my_rank) {  // on the same process
      MeshBlock *neighbor = pmb->pmy_mesh->FindMeshBlock(nb.snb.gid);
      auto exchanger = static_cast<T *>(
          neighbor->pimpl->GetExchanger(MessageTraits<T>::name));
      exchanger->recv_buffer_[nb.targetid].swap(Base::send_buffer_[nb.bufid]);
      exchanger->SetBoundaryStatus(nb.targetid, BoundaryStatus::arrived);
    }
#ifdef MPI_PARALLEL
    else {  // MPI
      int tag = create_mpi_tag(nb.snb.lid, nb.targetid, MessageTraits<T>::name);
      int ssize = Base::send_buffer_[nb.bufid].size();
      MPI_Isend(Base::send_buffer_[nb.bufid].data(), ssize,
                MessageTraits<T>::mpi_type, nb.snb.rank, tag, mpi_comm_,
                &req_mpi_send_[nb.bufid]);
    }
#endif
  }

  int rsize, tag;

#ifdef MPI_PARALLEL
  MPI_Status status;
  for (auto &nb : pmb->pbval->neighbor) {
    if (nb.snb.rank == Globals::my_rank) continue;  // local boundary received

    int tag = create_mpi_tag(pmb->lid, nb.bufid, Messeger<T>::name);

    MPI_Probe(nb.snb.rank, tag, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, Messeger<T>::mpi_type, &rsize);

    recv_buffer_[nb.bufid].resize(rsize);
    MPI_Irecv(recv_buffer_[nb.bufid].data(), rsize, MessegeTraits<T>::mpi_type,
              nb.snb.rank, tag, mpi_comm_, &req_mpi_recv_[nb.bufid]);
  }
#endif
}

template <typename T>
void NeighborExchanger<T>::ClearBuffer(MeshBlock const *pmb) {
  for (auto &nb : pmb->pbval->neighbor) {
#ifdef MPI_PARALLEL
    if (nb.snb.rank != Globals::my_rank)
      MPI_Wait(&req_mpi_send_[nb.bufid], MPI_STATUS_IGNORE);
#endif
  }
}

#endif  // SRC_EXCHANGER_NEIGHBOR_EXCHANGER_HPP_
