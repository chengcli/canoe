// C/C++
#include <algorithm>

// canoe
#include <configure.hpp>

// exchanger
#include <exchanger/linear_exchanger.hpp>

// harp
#include "radiation_band.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>

MPI_Datatype MessageTraits<RadiationBand>::mpi_type = MPI_ATHENA_REAL;

#endif  // MPI_PARALLEL

void RadiationBand::PackTemperature() {
  send_buffer_[0].resize(temf_.size());
  send_buffer_[0].swap(temf);
}

bool RadiationBand::UnpackTemperature() {
  int nblocks = 1;
#ifdef MPI_PARALLEL
  MPI_Comm_size(mpi_comm_, &nblocks);
#endif  // MPI_PARALLEL

#ifdef RT_DISORT
  auto& ds = psolver->ds_;

  for (int i = 0; i < (iu - il + 1) * nblocks; ++i) {
    int m = i / (iu - il + 1);
    ds.temper[m * (iu - il) + i % (iu - il + 1)] = recv_buffer_[i];
  }

  std::reverse(ds.temper, ds.temper + ds.nlyr + 1);
#endif  // RT_DISORT

  return true;
}

void RadiationBand::PackSpectralGrid(int b) {
  send_buffer_[1].resize(nlayer * (nphase_moments_ + 3));

  auto buf = send_buffer_[1].data();
  int nlayer = GetNumLayers();

  for (int i = 0; i < nlayer; ++i) {
    int i1 = i + NGHOST;
    *(buf++) = tau_(b, i1);
    *(buf++) = ssa_(b, i1);
    for (int j = 0; j <= nphase_moments_; ++j) {
      *(buf++) = pmom_(b, i1, j);
    }
  }
}

bool RadiationBand::UnpackSpectralGrid() {
  int nblocks = 1;
  int nlayer = GetNumLayers();
  int npmom = nphase_moments_;

#ifdef MPI_PARALLEL
  MPI_Comm_size(mpi_comm_, &nblocks);
#endif  // MPI_PARALLEL

  auto buf = recv_buffer_[1].data();

#ifdef RT_DISORT
  auto& ds = me->psolver->ds_;

  for (int n = 0; n < nblocks; ++n) {
    for (int i = 0; i < nlayer; ++i) {
      ds.dtauc[n * nlayer + i] = *(buf++);
      ds.ssalb[n * nlayer + i] = *(buf++);
      for (int j = 0; j <= npmom; ++j)
        ds.pmom[n * nlayer * ds.nmom_nstr + i * ds.nmom_nstr + j] = *(buf++);
      for (int j = npmom + 1; j < ds.nmom_nstr; ++j)
        ds.pmom[n * nlayer * ds.nmom_nstr + i * ds.nmom_nstr + j] = 0.;
    }
  }

  // absorption
  std::reverse(ds.dtauc, ds.dtauc + ds.nlyr);

  // single scatering albedo
  std::reverse(ds.ssalb, ds.ssalb + ds.nlyr);

  // phase function moments
  std::reverse(ds.pmom, ds.pmom + ds.nlyr * (ds.nmom_nstr + 1));
  for (int i = 0; i < ds.nlyr; ++i) {
    std::reverse(ds.pmom + i * (ds.nmom_nstr + 1),
                 ds.pmom + (i + 1) * (ds.nmom_nstr + 1));
  }
#endif  // RT_DISORT

  return true;
}

//! \bug only work for one block per process
void RadiationBand::Transfer(MeshBlock const* pmb, int n) {
  int nblocks = 1;
  int nlayer = GetNumLayers();
  int size = send_buffer_[n].size();

#ifdef MPI_PARALLEL
  MPI_Comm_size(mpi_comm_, &nblocks);
#endif  // MPI_PARALLEL

  if (n == 0) {
    recv_buffer_[0].resize(nblocks * size);
  } else if (n == 1) {
    recv_buffer_[1].resize(nblocks * nlayer * (nphase_moments_ + 3));
  }

#ifdef MPI_PARALLEL
  MPI_Allgather(&send_buffer_[n], size, MessageTraits<RadiationBand>::mpi_type,
                &recv_buffer_[n], size, MessageTraits<RadiationBand>::mpi_type,
                mpi_comm_);
#else
  recv_buffer_[n].swap(send_buffer_[n]);
#endif  // MPI_PARALLEL
}
