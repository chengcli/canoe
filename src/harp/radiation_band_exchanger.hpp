// canoe
#include <configure.hpp>

// exchanger
#include <exchanger/column_exchanger.hpp>

// harp
#include "radiation_band.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// Specialization for ParticleData Exchanger
template <>
struct Message<RadiationBand> {
  constexpr static int num_buffers = 1;
  constexpr static std::string name = "RadiationBand";

#ifdef MPI_PARALLEL
  static MPI_Datatype mpi_type;
#endif  // MPI_PARALLEL
};

template <>
ColumnExchanger<RadiationBand>::PackData() {
  T* buf = send_buffer_[0];
  int npmom = pmy_->GetNumPhaseMoments();

  for (int i = 0; i < nlayer; ++i) {
    *(buf++) = pmy_->tau_[i];
    *(buf++) = pmy_->ssa_[i];
    for (int j = 0; j < npmom; ++j) {
      *(buf++) = pmom[i * npmom + j];
    }
  }
}

template <>
ColumnExchanger<RadiationBand>::UnpackData() {
  npmom_max = std::max(npmom_max, npmom);
  T* buf = recv_buffer_[0];

  for (int n = 0; n < nblocks; ++n) {
    for (int i = 0; i < nlayer; ++i) {
      tau[n * nlayer + i] = *(buf++);
      ssa[n * nlayer + i] = *(buf++);
      for (int j = 0; j < npmom; ++j)
        pmom[n * nlayer * npmom_max + i * npmom_max + j] = *(buf++);
      for (int j = npmom; j < npmom_max; ++j)
        pmom[n * nlayer * npmom_max + i * npmom_max + j] = 0.;
    }
  }
}

template <>
ColumnExchanger<RadiationBand>::ClearBoundary() {
#ifdef MPI_PARALLEL
  MPI_Allgather(send_buffer_[0], send_size_, MPI_ATHENA_REAL, recv_buffer_[0],
                send_size_, MPI_ATHENA_REAL, mpi_comm_);
#else
  memcpy(recv_buffer_[0], send_buffer_[0], send_size_ * sizeof(Real));
#endif
}
