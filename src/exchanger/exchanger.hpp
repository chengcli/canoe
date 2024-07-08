#ifndef SRC_EXCHANGER_EXCHANGER_HPP_
#define SRC_EXCHANGER_EXCHANGER_HPP_

// C/C++
#include <array>
#include <vector>

// Athena++
#include <athena/athena.hpp>
#include <athena/bvals/bvals.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <configure.hpp>
#include <impl.hpp>
#include <virtual_groups.hpp>

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif  // MPI_PARALLEL

namespace ExchangeUtils {

int create_mpi_tag(int lid, int tid, std::string name);

//! find bottom neighbor block
NeighborBlock const *find_bot_neighbor(MeshBlock const *pmb);

//! find top neighbor block
NeighborBlock const *find_top_neighbor(MeshBlock const *pmb);

//! find left neighbor block
NeighborBlock const *find_left_neighbor(MeshBlock const *pmb);

//! find right neighbor block
NeighborBlock const *find_right_neighbor(MeshBlock const *pmb);

//! find back neighbor block
NeighborBlock const *find_back_neighbor(MeshBlock const *pmb);

//! find front neighbor block
NeighborBlock const *find_front_neighbor(MeshBlock const *pmb);

//! find neighbors in one coordinate direction
void find_neighbors(MeshBlock const *pmb, CoordinateDirection dir,
                    NeighborBlock *bblock, NeighborBlock *tblock);

}  // namespace ExchangeUtils

#ifndef MPI_PARALLEL

using MPI_Op = int;
extern int MPI_SUM;

#endif  // NOT MPI_PARALLEL

class ExchangerBase : public NamedGroup {
 public:  // constructor and destructor
  ExchangerBase(std::string name);
  virtual ~ExchangerBase();

  virtual void Regroup(MeshBlock const *pmb, CoordinateDirection dir) {}
  int GetRankInGroup() const;
  int GetGroupSize() const;
  void SetBlockID(int id) { gid_ = id; }

  std::vector<SimpleNeighborBlock> neighbor;

 protected:
  void setColor(MeshBlock const *pmb, CoordinateDirection dir);

#ifdef MPI_PARALLEL
  MPI_Comm mpi_comm_;
#endif

  //! \brief global block id associated with this exchanger
  int gid_;

  //! \brief MPI color of each block (nrank)
  std::vector<int> color_;

  //! \brief MPI rank of the bottom of each block (nrank)
  std::vector<int> brank_;
};

template <typename T, int N>
class Exchanger : public ExchangerBase {
 public:  // constructor and destructor
  using BufferType = std::vector<T>;

  std::array<BufferType, N> send_buffer;
  std::array<BufferType, N> recv_buffer;

#ifdef MPI_PARALLEL
  static MPI_Datatype mpi_type;
  using ExchangerBase::mpi_comm_;
#endif

  Exchanger(std::string name) : ExchangerBase(name) {
    for (int n = 0; n < N; ++n) {
#ifdef MPI_PARALLEL
      req_mpi_send_[n] = MPI_REQUEST_NULL;
      req_mpi_recv_[n] = MPI_REQUEST_NULL;
#endif  // MPI_PARALLEL
      status_flag_[n] = BoundaryStatus::waiting;
    }
  }

  virtual ~Exchanger() {
    for (int n = 0; n < N; ++n) {
#ifdef MPI_PARALLEL
      req_mpi_send_[n] = MPI_REQUEST_NULL;
      req_mpi_recv_[n] = MPI_REQUEST_NULL;
#endif  // MPI_PARALLEL
      status_flag_[n] = BoundaryStatus::waiting;
    }
  }

  //! \brief Send and receive data
  //! \bug only work for one block per process
  virtual void GatherAll(MeshBlock const *pmb) {
    int nblocks = GetGroupSize();

    if (nblocks > 1) {
#ifdef MPI_PARALLEL
      for (int n = 0; n < N; ++n) {
        int size = send_buffer[n].size();
        MPI_Allgather(send_buffer[n].data(), size, mpi_type,
                      recv_buffer[n].data(), size, mpi_type, mpi_comm_);
      }
#endif  // MPI_PARALLEL
    } else {
      for (int n = 0; n < N; ++n) recv_buffer[n].swap(send_buffer[n]);
    }
  }

  //! \brief Send and receive data
  virtual void ReduceAll(MeshBlock const *pmb, MPI_Op op) {
#ifdef MPI_PARALLEL
    // sum over all ranks
    for (int n = 0; n < N; ++n) {
      MPI_Allreduce(send_buffer[n].data(), recv_buffer[n].data(),
                    send_buffer[n].size(), mpi_type, op, mpi_comm_);
    }
#endif
  }

  //! \brief Point-to-Point communication
  virtual void RecvFrom(SimpleNeighborBlock const &snb) {
    if (snb.rank != Globals::my_rank) {  // MPI boundary
#ifdef MPI_PARALLEL
      MPI_Status status;
      for (int n = 0; n < N; ++n) {
        int tag = ExchangeUtils::create_mpi_tag(gid_, snb.gid,
                                                GetName() + std::to_string(n));
        MPI_Recv(recv_buffer[n].data(), recv_buffer[n].size(), mpi_type,
                 snb.rank, tag, mpi_comm_, &status);
        SetBoundaryStatus(n, BoundaryStatus::arrived);
      }
#endif
    }  // nothing to do for local boundary (taken care of in SendTo)
  }

  //! \brief Point-to-Point communication
  virtual void SendTo(MeshBlock const *pmb, SimpleNeighborBlock const &snb) {
    if (snb.rank != Globals::my_rank) {  // MPI boundary
#ifdef MPI_PARALLEL
      for (int n = 0; n < N; ++n) {
        int tag = ExchangeUtils::create_mpi_tag(snb.gid, gid_,
                                                GetName() + std::to_string(n));
        MPI_Isend(send_buffer[n].data(), send_buffer[n].size(), mpi_type,
                  snb.rank, tag, mpi_comm_, &req_mpi_send_[n]);
      }
#endif
    } else {  // local boundary
      MeshBlock *pbl = pmb->pmy_mesh->FindMeshBlock(snb.gid);
      auto pexv = static_cast<Exchanger<T, N> *>(
          pbl->pimpl->FindExchanger(GetName().c_str()));
      for (int n = 0; n < N; ++n) {
        pexv->recv_buffer[n].swap(send_buffer[n]);
        pexv->SetBoundaryStatus(n, BoundaryStatus::arrived);
      }
    }
  }

  void ClearStatus() {
    for (int n = 0; n < N; ++n) {
      status_flag_[n] = BoundaryStatus::waiting;
    }
  }

  //! \brief Clear buffer
  void ClearBuffer() {
#ifdef MPI_PARALLEL
    for (int n = 0; n < N; ++n) {
      MPI_Wait(&req_mpi_send_[n], MPI_STATUS_IGNORE);
    }
#endif  // MPI_PARALLEL
  }

  //! \brief Set the boundary status
  void SetBoundaryStatus(int bid, BoundaryStatus status) {
    status_flag_[bid] = status;
  }

  BoundaryStatus GetBoundaryStatus(int bid) const { return status_flag_[bid]; }

 protected:
  std::array<enum BoundaryStatus, N> status_flag_;

#ifdef MPI_PARALLEL
  std::array<MPI_Request, N> req_mpi_send_;
  std::array<MPI_Request, N> req_mpi_recv_;
#endif  // MPI_PARALLEL
};

#ifdef MPI_PARALLEL
template <typename T, int N>
MPI_Datatype Exchanger<T, N>::mpi_type = MPI_ATHENA_REAL;
#endif  // MPI_PARALLEL

template <typename T, int N>
class LinearExchanger : public Exchanger<T, N> {
 public:  // constructor and destructor
  using Base = Exchanger<T, N>;

#ifdef MPI_PARALLEL
  using Base::mpi_comm_;
#endif  // MPI_PARALLEL

  LinearExchanger(std::string name) : Base(name) {}

  void Regroup(MeshBlock const *pmb, CoordinateDirection dir) override {
    std::fill(Base::color_.begin(), Base::color_.end(), -1);
    Base::setColor(pmb, dir);
  }
};

template <typename T, int N>
class PlanarExchanger : public Exchanger<T, N> {
 public:  // constructor and destructor
  using Base = Exchanger<T, N>;

#ifdef MPI_PARALLEL
  using Base::mpi_comm_;
#endif  // MPI_PARALLEL

  PlanarExchanger(std::string name) : Base(name) {}

  void Regroup(MeshBlock const *pmb, CoordinateDirection dir) override {
    std::fill(Base::color_.begin(), Base::color_.end(), -1);

    if (dir == X1DIR) {
      Base::setColor(pmb, X2DIR);
      Base::setColor(pmb, X3DIR);
    } else if (dir == X2DIR) {
      Base::setColor(pmb, X1DIR);
      Base::setColor(pmb, X3DIR);
    } else {
      Base::setColor(pmb, X1DIR);
      Base::setColor(pmb, X2DIR);
    }
  }
};

#endif  // SRC_EXCHANGER_EXCHANGER_HPP_
