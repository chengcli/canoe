// C/C++
#include <algorithm>

// canoe
#include <configure.hpp>

// exchanger
#include <exchanger/exchanger.hpp>

// harp
#include "radiation_band.hpp"
#include "rt_solvers.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif  // MPI_PARALLEL

void RadiationBand::PackTemperature() {
  send_buffer_[0].resize(temf_.size());
  send_buffer_[0].swap(temf_);
}

bool RadiationBand::UnpackTemperature(void *arg) {
  int nblocks = 1;
  int nlevels = temf_.size();
  int nlayers = GetNumLayers();

#ifdef MPI_PARALLEL
  int is_initialized;
  MPI_Initialized(&is_initialized);

  if (is_initialized) MPI_Comm_size(mpi_comm_, &nblocks);
#endif  // MPI_PARALLEL

#ifdef RT_DISORT
  disort_state *ds = static_cast<disort_state *>(arg);

  for (int n = 0; n < nblocks; ++n) {
    for (int i = 0; i <= nlayers; ++i) {
      ds->temper[n * nlayers + i] = recv_buffer_[0][n * nlevels + i + NGHOST];
    }
  }

  std::reverse(ds->temper, ds->temper + ds->nlyr + 1);
#endif  // RT_DISORT

  return true;
}

void RadiationBand::PackSpectralGrid(int b) {
  int nlayers = GetNumLayers();
  int npmom = GetNumPhaseMoments();
  send_buffer_[1].resize(nlayers * (npmom + 3));

  auto buf = send_buffer_[1].data();

  for (int i = 0; i < nlayers; ++i) {
    int i1 = i + NGHOST;
    *(buf++) = tau_(b, i1);
    *(buf++) = ssa_(b, i1);
    for (int j = 0; j <= npmom; ++j) {
      *(buf++) = pmom_(b, i1, j);
    }
  }
}

bool RadiationBand::UnpackSpectralGrid(void *arg) {
  int nblocks = 1;
  int nlayers = GetNumLayers();
  int npmom = GetNumPhaseMoments();

#ifdef MPI_PARALLEL
  int is_initialized;
  MPI_Initialized(&is_initialized);

  if (is_initialized) MPI_Comm_size(mpi_comm_, &nblocks);
#endif  // MPI_PARALLEL

  auto buf = recv_buffer_[1].data();

#ifdef RT_DISORT
  disort_state *ds = static_cast<disort_state *>(arg);

  for (int n = 0; n < nblocks; ++n) {
    for (int i = 0; i < nlayers; ++i) {
      ds->dtauc[n * nlayers + i] = *(buf++);
      ds->ssalb[n * nlayers + i] = *(buf++);
      for (int j = 0; j <= npmom; ++j) {
        ds->pmom[n * nlayers * (ds->nmom_nstr + 1) + i * (ds->nmom_nstr + 1) +
                 j] = *(buf++);
      }
      for (int j = npmom + 1; j <= ds->nmom_nstr; ++j)
        ds->pmom[n * nlayers * (ds->nmom_nstr + 1) + i * (ds->nmom_nstr + 1) +
                 j] = 0.;
    }
  }

  // absorption
  std::reverse(ds->dtauc, ds->dtauc + ds->nlyr);

  // single scatering albedo
  std::reverse(ds->ssalb, ds->ssalb + ds->nlyr);

  // phase function moments
  std::reverse(ds->pmom, ds->pmom + ds->nlyr * (ds->nmom_nstr + 1));
  for (int i = 0; i < ds->nlyr; ++i) {
    std::reverse(ds->pmom + i * (ds->nmom_nstr + 1),
                 ds->pmom + (i + 1) * (ds->nmom_nstr + 1));
  }
#endif  // RT_DISORT

  return true;
}

//! \bug only work for one block per process
void RadiationBand::Transfer(MeshBlock const *pmb, int n) {
  int nblocks = 1;
  int nlayers = GetNumLayers();
  int npmom = GetNumPhaseMoments();
  int size = send_buffer_[n].size();

#ifdef MPI_PARALLEL
  int is_initialized;
  MPI_Initialized(&is_initialized);

  if (is_initialized) MPI_Comm_size(mpi_comm_, &nblocks);
#endif  // MPI_PARALLEL

  if (n == 0) {
    recv_buffer_[0].resize(nblocks * size);
  } else if (n == 1) {
    recv_buffer_[1].resize(nblocks * nlayers * (npmom + 3));
  }

  if (nblocks > 1) {
#ifdef MPI_PARALLEL
    MPI_Allgather(&send_buffer_[n], size, MPI_ATHENA_REAL, &recv_buffer_[n],
                  size, MPI_ATHENA_REAL, mpi_comm_);
#endif  // MPI_PARALLEL
  } else {
    recv_buffer_[n].swap(send_buffer_[n]);
  }
}
