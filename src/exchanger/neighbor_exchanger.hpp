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
void NeighborExchanger<T>::SendBuffer() const {
  auto pmb = getMeshBlock();

  for (auto &nb : pmb->pbval->neighbor) {
    if (nb.snb.rank == Globals::my_rank) {  // on the same process
      MeshBlock *neighbor = pmb->pmy_mesh->FindMeshBlock(nb.snb.gid);
      auto &exchanger = static_cast<NeighborExchanger<T> *>(
          neighbor->pimpl->GetExchanger(Message<T>::name));
      exchanger->SetRecvBuffer(nb.targetid, send_buffer_[nb.bufid]);
      exchanger->SetBoundaryStatus(nb.targetid, BoundaryStatus::arrived);
    }
#ifdef MPI_PARALLEL
    else {  // MPI
      int tag = create_mpi_tag(nb.snb.lid, nb.targetid, Message<T>::name);
      int ssize = send_buffer_[nb.bufid].size();
      MPI_Isend(send_buffer_[nb.bufid], ssize, Message<T>::mpi_type,
                nb.snb.rank, tag, mpi_comm_, &req_mpi_send_[nb.bufid]);
    }
#endif
  }
}

template <typename T>
void NeighborExchanger<T>::RecvRuffer() {
  int rsize, tag;
  auto pmb = getMeshBlock();

#ifdef MPI_PARALLEL
  MPI_Status status;
  for (auto &nb : pmb->pbval->neighbor) {
    if (nb.snb.rank == Globals::my_rank) continue;  // local boundary received

    if (status_flag_[nb.bufid] == BoundaryStatus::waiting) {
      int tag = create_mpi_tag(pmb->lid, nb.bufid, Messeger<T>::name);

      MPI_Probe(nb.snb.rank, tag, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, Messeger<T>::mpi_type, &rsize);

      recv_buffer_[nb.bufid].resize(rsize);
      MPI_Irecv(recv_buffer_[nb.bufid], rsize, Messeger<T>::mpi_type,
                nb.snb.rank, tag, mpi_comm_, &req_mpi_recv_[nb.bufid]);
    }
  }
#endif
}

template <typename T>
void NeighborExchanger<T>::ClearBoundary() {
  auto pmb = getMeshBlock();

  for (auto &nb : pmb->pbval->neighbor) {
    status_flag_[nb.bufid] = BoundaryStatus::waiting;

#ifdef MPI_PARALLEL
    if (nb.snb.rank != Globals::my_rank)
      MPI_Wait(&req_mpi_send_[nb.bufid], MPI_STATUS_IGNORE);
#endif
  }
}

#endif  // SRC_EXCHANGER_NEIGHBOR_EXCHANGER_HPP_
