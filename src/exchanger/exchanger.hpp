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

class ExchangerBase {
 public:
  ExchangerBase() {
#ifdef MPI_PARALLEL
    mpi_comm_ = MPI_COMM_WORLD;
#endif
  }

  virtual ExchangerBase() {}

  virtual void SendBuffer() const {}
  virtual void RecvBuffer() {}

  virtual void ClearBuffer() = 0;

  virtual void PackData() = 0;
  virtual bool UnpackData() = 0;

 protected:
#ifdef MPI_PARALLEL
  MPI_Comm mpi_comm_;
#endif
};

template <typename T>
class Exchanger : public ExchangerBase {
 public:  // constructor and destructor
  using DataType = MessageTraits<T>::DataType;
  using BufferType = std::vector<DataType>;

  Exchanger(T *me) : phost_(me) {
    for (int n = 0; n < MessageTraits<T>::num_buffers; ++n) {
      status_flag_[n] = BoundaryStatus::waiting;

#ifdef MPI_PARALLEL
      req_mpi_send_[n] = MPI_REQUEST_NULL;
      req_mpi_recv_[n] = MPI_REQUEST_NULL;
#endif
    }
  }

  virtual ~Exchanger() {
    for (int n = 0; n < MessageTraits<T>::num_buffers; ++n) {
      status_flag_[n] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
      req_mpi_send_[n] = MPI_REQUEST_NULL;
      req_mpi_recv_[n] = MPI_REQUEST_NULL;

      if (mpi_comm_ != MPI_COMM_WORLD) {
        MPI_Comm_free(&mpi_comm_);
      }
#endif
    }
  }

  void SetBoundaryStatus(int bid, BoundaryStatus status) {
    status_flag_[bid] = status;
  }

  T const *Me() const { return phost_; }

 protected:
  enum BoundaryStatus status_flag_[MessageTraits<T>::num_buffers];

  BufferType send_buffer_[MessageTraits<T>::num_buffers];
  BufferType recv_buffer_[MessageTraits<T>::num_buffers];

#ifdef MPI_PARALLEL
  MPI_Request req_mpi_send_[MessageTraits<T>::num_buffers];
  MPI_Request req_mpi_recv_[MessageTraits<T>::num_buffers];
#endif

 private:
  T const *phost_;
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
void find_neighbors(MeshBlock const *pmb, CoordianteID dir,
                    NeighborBlock *bblock, NeighborBlock *tblock);
}  // namespace ExchangerHelper

template <typename T>
class NeighborExchanger : public Exchanger<T> {
 public:
  NeighborExchanger(T *me);
  virtual void SendBuffer() const override;
  virtual void RecvRuffer() override;
  virtual void ClearBoundary() override;
};

template <typename T>
class LinearExchanger : public Exchanger<T> {
 public:  // constructor and destructor
  LinearExchanger(T *me);

 public:  // member functions
  int GetRankInGroup() const;
  void Regroup(CoordinateID dir);

 public:  // exchanger functions
  virtual void ClearBoundary() override;
  virtual void PackData() override;
  virtual bool UnpackData() override;

 protected:
  //! \brief MPI color of each block
  std::vector<int> color_;

  //! \brief MPI rank of the bottom of each block
  std::vector<int> brank_;
};

template <typename T>
class PlanarExchanger : public Exchanger<T> {
 public:
  PlanarExchanger(T *me);
};

template <typename T>
class VolumetricExchanger : public Exchanger<T> {
 public:
  VolumetricExchanger(T *me);
};

#endif  // SRC_EXCHANGER_EXCHANGER_HPP_
