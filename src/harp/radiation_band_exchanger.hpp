// canoe
#include <configure.hpp>

// exchanger
#include <exchanger/linear_exchanger.hpp>

// harp
#include "radiation_band.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

template <>
void LinearExchanger<RadiationBand>::PackData(void* args) {
  int nlayer = Me()->GetNumLayers();
  int npmom = Me()->GetNumPhaseMoments();

  // first send buffer is temperature at levels
  send_buffer_[0].swap(temf);

  auto& tau = Me()->tau_;
  auto& ssa = Me()->ssa_;
  auto& pmom = Me()->pmom_;

  // second send buffer is tau, ssa, and pmom at layers
  auto buf = send_buffer_[1].data();

  for (int i = 0; i < nlayer; ++i) {
    *(buf++) = tau(n, i);
    *(buf++) = ssa(n, i);
    for (int j = 0; j < npmom; ++j) {
      *(buf++) = pmom[i * npmom + j];
    }
  }
}

template <>
bool LinearExchanger<RadiationBand>::UnpackData() {
  npmom_max = std::max(npmom_max, npmom);

  auto buf = recv_buffer_[0].data();

#ifdef DISORT
  auto ds = Me()->psolver->ds_;

  for (int n = 0; n < nblocks; ++n) {
    for (int i = 0; i < nlayer; ++i) {
      ds->dtauc[n * nlayer + i] = *(buf++);
      ds->ssalb[n * nlayer + i] = *(buf++);
      for (int j = 0; j < npmom; ++j)
        ds->pmom[n * nlayer * npmom_max + i * npmom_max + j] = *(buf++);
      for (int j = npmom; j < npmom_max; ++j)
        ds->pmom[n * nlayer * npmom_max + i * npmom_max + j] = 0.;
    }
  }
#endif

  return true;
}

template <>
void LinearExchanger<RadiationBand>::Execute() {
  for (int n = 0; n < MessageTraits<RadiationBand>::num_buffers; ++n) {
    int size = send_buffer_[n].size();

#ifdef MPI_PARALLEL
    MPI_Allgather(&send_buffer_[n], size,
                  MessageTraits<RadiationBand>::mpi_type, &recv_buffer_[n],
                  size, MessageTraits<RadiationBand>::mpi_trype, mpi_comm_);
#else
    memcpy(&recv_buffer_[n], &send_buffer_[n], size * sizeof(Real));
#endif
  }
}
