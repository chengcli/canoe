#ifndef SRC_EXCHANGER_EXCHANGER_HPP_
#define SRC_EXCHANGER_EXCHANGER_HPP_

// Athena++
#include <athena/bvals/bvals.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <configure.hpp>
#include <impl.hpp>

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif  // MPI_PARALLEL

class NeighborBlock;
class MeshBlock;

class ExchangerBase {
 public:  // constructor and destructor
  ExchangerBase();
  virtual ~ExchangerBase() {}

  //! \brief Pack data to send buffer
  virtual void PackData(MeshBlock const *pmb) {}

  //! \brief Unpack data from receive buffer
  virtual bool UnpackData(MeshBlock const *pmb) { return true; }

  //! \brief Send and receive data
  virtual void Transfer(MeshBlock const *pmb, int n = -1) {}

  virtual void Regroup(MeshBlock const *pmb, CoordinateDirection dir) {}

  int GetRankInGroup() const;

 protected:
  void setColor(MeshBlock const *pmb, CoordinateDirection dir);

#ifdef MPI_PARALLEL
  MPI_Comm mpi_comm_;
#endif

  //! \brief MPI color of each block
  std::vector<int> color_;

  //! \brief MPI rank of the bottom of each block
  std::vector<int> brank_;
};

template <typename T, int N>
class Exchanger : public ExchangerBase {
 public:  // constructor and destructor
  using BufferType = std::vector<T>;

#ifdef MPI_PARALLEL
  using ExchangerBase::mpi_comm_;
#endif

  Exchanger() {
#ifdef MPI_PARALLEL
    for (int n = 0; n < N; ++n) {
      req_mpi_send_[n] = MPI_REQUEST_NULL;
      req_mpi_recv_[n] = MPI_REQUEST_NULL;
    }
#endif  // MPI_PARALLEL
  }

  virtual ~Exchanger() {
#ifdef MPI_PARALLEL
    for (int n = 0; n < N; ++n) {
      req_mpi_send_[n] = MPI_REQUEST_NULL;
      req_mpi_recv_[n] = MPI_REQUEST_NULL;
    }

    if (mpi_comm_ != MPI_COMM_WORLD) {
      MPI_Comm_free(&mpi_comm_);
    }
#endif  // MPI_PARALLEL
  }

  //! \brief Clear buffer
  virtual void ClearBuffer(MeshBlock const *pmb) {
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

 protected:
  std::array<enum BoundaryStatus, N> status_flag_;
  std::array<BufferType, N> send_buffer_;
  std::array<BufferType, N> recv_buffer_;

#ifdef MPI_PARALLEL
  std::array<MPI_Request, N> req_mpi_send_;
  std::array<MPI_Request, N> req_mpi_recv_;
#endif  // MPI_PARALLEL
};

namespace ExchangerHelper {

extern int mpi_tag_ub;
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

}  // namespace ExchangerHelper

template <typename T, int N>
class LinearExchanger : public Exchanger<T, N> {
 public:  // constructor and destructor
  using Base = Exchanger<T, N>;

#ifdef MPI_PARALLEL
  using Base::mpi_comm_;
#endif  // MPI_PARALLEL

  LinearExchanger() {}

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

  PlanarExchanger() {}

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
