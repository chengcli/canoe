#ifndef SRC_EXCHANGER_EXCHANGER_HPP_
#define SRC_EXCHANGER_EXCHANGER_HPP_

// Athena++
#include <athena/mesh/mesh.hpp>

// canoe
#include <configure.hpp>
#include <impl.hpp>

// exchanger
#include "message_traits.hpp"

class NeighborBlock;
class MeshBlock;

class ExchangerBase {
 public:  // constructor and destructor
  ExchangerBase() {
#ifdef MPI_PARALLEL
    mpi_comm_ = MPI_COMM_WORLD;
#endif
  }
  virtual ~ExchangerBase() {}

  //! \brief Pack data to send buffer
  virtual void PackData(MeshBlock const *pmb) {}

  //! \brief Unpack data from receive buffer
  virtual bool UnpackData(MeshBlock const *pmb) { return true; }

  //! \brief Send and receive data
  virtual void Transfer(MeshBlock const *pmb, int n = -1) = 0;

 protected:
#ifdef MPI_PARALLEL
  MPI_Comm mpi_comm_;
#endif
};

template <typename T>
class Exchanger : public ExchangerBase {
 public:  // constructor and destructor
  using DataType = typename MessageTraits<T>::DataType;
  using BufferType = std::vector<DataType>;

  Exchanger() {
    for (int n = 0; n < MessageTraits<T>::num_buffers; ++n) {
#ifdef MPI_PARALLEL
      req_mpi_send_[n] = MPI_REQUEST_NULL;
      req_mpi_recv_[n] = MPI_REQUEST_NULL;
#endif  // MPI_PARALLEL
    }
  }

  virtual ~Exchanger() {
    for (int n = 0; n < MessageTraits<T>::num_buffers; ++n) {
#ifdef MPI_PARALLEL
      req_mpi_send_[n] = MPI_REQUEST_NULL;
      req_mpi_recv_[n] = MPI_REQUEST_NULL;

      if (mpi_comm_ != MPI_COMM_WORLD) {
        MPI_Comm_free(&mpi_comm_);
      }
#endif  // MPI_PARALLEL
    }
  }

  //! \brief Clear buffer
  virtual void ClearBuffer(MeshBlock const *pmb) {
    for (int n = 0; n < MessageTraits<T>::num_buffers; ++n) {
#ifdef MPI_PARALLEL
      MPI_Wait(&req_mpi_send_[n], MPI_STATUS_IGNORE);
#endif  // MPI_PARALLEL
    }
  }

  //! \brief Set the boundary status
  void SetBoundaryStatus(int bid, BoundaryStatus status) {
    status_flag_[bid] = status;
  }

 protected:
  enum BoundaryStatus status_flag_[MessageTraits<T>::num_buffers];
  BufferType send_buffer_[MessageTraits<T>::num_buffers];
  BufferType recv_buffer_[MessageTraits<T>::num_buffers];

#ifdef MPI_PARALLEL
  MPI_Request req_mpi_send_[MessageTraits<T>::num_buffers];
  MPI_Request req_mpi_recv_[MessageTraits<T>::num_buffers];
#endif  // MPI_PARALLEL
};

namespace ExchangerHelper {

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

template <typename T>
class NeighborExchanger : public Exchanger<T> {
 public:
  using Base = Exchanger<T>;
  using Base::recv_buffer_;
  using Base::req_mpi_recv_;
  using Base::req_mpi_send_;
  using Base::send_buffer_;
  using ExchangerBase::mpi_comm_;

  NeighborExchanger() {}

  //! \brief Send and receive data
  virtual void Transfer(MeshBlock const *pmb, int n = -1) override;

  //! \brief Clear buffer
  virtual void ClearBuffer(MeshBlock const *pmb) override;
};

template <typename T>
class LinearExchanger : public Exchanger<T> {
 public:  // constructor and destructor
  using Base = Exchanger<T>;
  using ExchangerBase::mpi_comm_;

  LinearExchanger();

 public:  // member functions
  int GetRankInGroup() const;
  void Regroup(MeshBlock const *pmb, CoordinateDirection dir);

 private:
  //! \brief MPI color of each block
  std::vector<int> color_;

  //! \brief MPI rank of the bottom of each block
  std::vector<int> brank_;
};

template <typename T>
class PlanarExchanger : public Exchanger<T> {
 public:
  PlanarExchanger();
};

#endif  // SRC_EXCHANGER_EXCHANGER_HPP_
