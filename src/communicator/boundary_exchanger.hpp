#ifndef SRC_COMMUNICATOR_BOUNDARY_EXCHANGER_HPP_
#define SRC_COMMUNICATOR_BOUNDARY_EXCHANGER_HPP_

// C/C++
#include <functional> // hash

// Athena++
#include <athena/bvals/bvals.hpp>
#include <athena/globals.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <configure.hpp>
#include <impl.hpp>

// communicator
#include "communicator.hpp"

namespace Globals {
extern int mpi_tag_ub;
}

// helper functions
inline int create_mpi_tag(int lid, int tid, std::string exchanger_name) {
  int tag = BoundaryBase::CreateBvalsMPITag(lid, tid, 0);

  std::string str = exchanger_name + std::to_string(tag);
  return std::hash<std::string>{}(str)%(Globals::mpi_tag_ub);
}

template <typename T>
class BoundaryExchanger {
 public:
  BoundaryExchanger() {
    for (int i = 0; i < Communicator<T>::max_neighbors; ++i) {
      status_flag_[i] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
      req_mpi_send_[i] = MPI_REQUEST_NULL;
      req_mpi_recv_[i] = MPI_REQUEST_NULL;
#endif
    }
  }

  virtual ~BoundaryExchanger() {}

  void SetRecvBuffer(int bid, std::vector<T> const &buffer) {
    recv_buffer_[bid] = buffer;
  }

  void SetBoundaryStatus(int bid, BoundaryStatus status) {
    status_flag_[bid] = status;
  }

  void SendBuffer() const {
    auto pmb = getMeshBlock();
    for (auto &nb : pmb->pbval->neighbor) {
      if (nb.snb.rank == Globals::my_rank) {  // on the same process
        MeshBlock *neighbor = pmb->pmy_mesh->FindMeshBlock(nb.snb.gid);
        auto &exchanger = static_cast<BoundaryExchanger<T> *>(
            neighbor->pimpl->GetExchanger(Communicator<T>::name));
        exchanger->SetRecvBuffer(nb.targetid, send_buffer_[nb.bufid]);
        exchanger->SetBoundaryStatus(nb.targetid, BoundaryStatus::arrived);
      }
#ifdef MPI_PARALLEL
      else {  // MPI
        int tag = create_mpi_tag(nb.snb.lid, nb.targetid, Communicator<T>::name);
        int ssize = send_buffer_[nb.bufid].size();
        MPI_Isend(send_buffer_[nb.bufid].data(), ssize,
                  Communicator<T>::mpi_type, nb.snb.rank, tag,
                  Communicator<T>::mpi_comm, &req_mpi_send_[nb.bufid]);
      }
#endif
    }
  }

  void RecvRuffer() {
    int rsize, tag;
    auto pmb = getMeshBlock();

#ifdef MPI_PARALLEL
    MPI_Status status;
    for (auto &nb : pmb->pbval->neighbor) {
      if (nb.snb.rank == Globals::my_rank) continue;  // local boundary received

      if (status_flag_[nb.bufid] == BoundaryStatus::waiting) {
        int tag = create_mpi_tag(pmb->lid, nb.bufid, Communicator<T>::name);

        MPI_Probe(nb.snb.rank, tag, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, Communicator<T>::mpi_type, &rsize);

        recv_buffer_[nb.bufid].resize(rsize);
        MPI_Irecv(recv_buffer_[nb.bufid].data(), rsize,
                  Communicator<T>::mpi_type, nb.snb.rank, tag,
                  Communicator<T>::mpi_comm, &req_mpi_recv_[nb.bufid]);
      }
    }
#endif
  }

  void ClearBoundary() {
    auto pmb = getMeshBlock();

    for (auto &nb : pmb->pbval->neighbor) {
      status_flag_[nb.bufid] = BoundaryStatus::waiting;
      send_buffer_[nb.bufid].clear();
      recv_buffer_[nb.bufid].clear();

#ifdef MPI_PARALLEL
      if (nb.snb.rank != Globals::my_rank)
        MPI_Wait(&req_mpi_send_[nb.bufid], MPI_STATUS_IGNORE);
#endif
    }
  }

  virtual void DetachTo(std::vector<T> &buffer) = 0;
  virtual bool AttachTo(typename Communicator<T>::container_t &container) = 0;

 protected:
  virtual MeshBlock const *getMeshBlock() const = 0;

  enum BoundaryStatus status_flag_[Communicator<T>::max_neighbors];
  std::vector<T> send_buffer_[Communicator<T>::max_neighbors];
  std::vector<T> recv_buffer_[Communicator<T>::max_neighbors];

#ifdef MPI_PARALLEL
  MPI_Request req_mpi_send_[Communicator<T>::max_neighbors];
  MPI_Request req_mpi_recv_[Communicator<T>::max_neighbors];
#endif
};

#endif  // SRC_COMMUNICATOR_BOUNDARY_EXCHANGER_HPP_
