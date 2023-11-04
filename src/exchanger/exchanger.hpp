#ifndef SRC_EXCHANGER_EXCHANGER_HPP_
#define SRC_EXCHANGER_EXCHANGER_HPP_

// C/C++
#include <functional>  // hash

// Athena++
#include <athena/bvals/bvals.hpp>
#include <athena/globals.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <configure.hpp>
#include <impl.hpp>

// exchanger
#include "communicator.hpp"

template <typename T>
class Exchanger {
 public:  // constructor and destructor
  using DataType = Message<T>::DataType;

  Exchanger(T *me) : pmy_(me) {
    for (int n = 0; n < Message<T>::num_buffers; ++n) {
      send_buffer_[n] = nullptr;
      recv_buffer_[n] = nullptr;
      status_flag_[n] = BoundaryStatus::waiting;

#ifdef MPI_PARALLEL
      req_mpi_send_[n] = MPI_REQUEST_NULL;
      req_mpi_recv_[n] = MPI_REQUEST_NULL;
#endif
    }
  }

  virtual ~Exchanger() {
    for (int n = 0; n < Message<T>::num_buffers; ++n) {
      if (send_buffer_[n] != nullptr) {
        delete[] send_buffer_[n];
      }

      if (recv_buffer_[n] != nullptr) {
        delete[] recv_buffer_[n];
      }

      status_flag_[i] = BoundaryStatus::waiting;
#ifdef MPI_PARALLEL
      req_mpi_send_[n] = MPI_REQUEST_NULL;
      req_mpi_recv_[n] = MPI_REQUEST_NULL;

      if (mpi_comm_ != MPI_COMM_WORLD) {
        MPI_Comm_free(&mpi_comm_);
      }
#endif
    }
  }

  void SetRecvBuffer(int bid, data_type *buffer) { recv_buffer_[bid] = buffer; }

  void SetBoundaryStatus(int bid, BoundaryStatus status) {
    status_flag_[bid] = status;
  }

  void ResizeSendBuffer(int id, int size) {
    if (send_buffer_[id] != nullptr) {
      delete[] send_buffer_[id];
    }

    send_buffer_[id] = new DataType[size];
  };

  void ResizeRecvBuffer(int id, int size) {
    if (recv_buffer_[id] != nullptr) {
      delete[] recv_buffer_[id];
    }

    recv_buffer_[id] = new DataType[size];
  };

  virtual void SendBuffer() const = 0;
  virtual void RecvRuffer() = 0;
  virtual void ClearBoundary() = 0;

  virtual void PackData() = 0;
  virtual bool UnpackData() = 0;

  // void reduceData23(Real *send, Real *recv);
  void GatherData(Real *send, Real *recv, int size) const;
  void GatherDataInPlace(Real *recv, int size) const;

 protected:
  T const *pmy_;

  enum BoundaryStatus status_flag_[Message<T>::num_buffers];
  DataType *send_buffer_[Message<T>::num_buffers];
  DataType *recv_buffer_[Message<T>::num_buffers];

#ifdef MPI_PARALLEL
  MPI_Request req_mpi_send_[Message<T>::num_buffers];
  MPI_Request req_mpi_recv_[Message<T>::num_buffers];
#endif

#ifdef MPI_PARALLEL
  MPI_Comm mpi_comm_;
#endif
};

template <typename T>
class NeighbourExchanger : public Exchanger<T> {
 public:
  NeighbourExchanger(T *me);
  virtual void SendBuffer() const override;
  virtual void RecvRuffer() override;
  virtual void ClearBoundary() override;
};

template <typename T>
class ColumnExchanger : public Exchanger<T> {
 public:
  ColumnExchanger(T *me) : Exchanger<T>(me) {}

  int GetRankInGroup() const;
  void Regroup(CoordinateID dir);

 protected:
  std::vector<int> color_;
  std::vector<int> brank_;
};

template <typename T>
class PlaneExchanger : public Exchanger<T> {
 public:
  PlaneExchanger(T *me);
};

#endif  // SRC_EXCHANGER_EXCHANGER_HPP_
