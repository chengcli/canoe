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

void RadiationBand::packTemperature() { pexv->send_buffer[0].swap(temf_); }

bool RadiationBand::unpackTemperature(void *arg) {
  int nblocks = pexv->GetGroupSize();
  int nlevels = temf_.size();
  int nlayers = GetNumLayers();

#ifdef RT_DISORT
  disort_state *ds = static_cast<disort_state *>(arg);

  for (int n = 0; n < nblocks; ++n) {
    for (int i = 0; i <= nlayers; ++i) {
      ds->temper[n * nlayers + i] =
          pexv->recv_buffer[0][n * nlevels + i + NGHOST];
    }
  }

  std::reverse(ds->temper, ds->temper + ds->nlyr + 1);
#endif  // RT_DISORT

  return true;
}

void RadiationBand::packSpectralProperties() {
  int nlayers = GetNumLayers();
  int npmom = GetNumPhaseMoments();

  auto buf = pexv->send_buffer[1].data();

  for (int b = 0; b < pgrid_->spec.size(); ++b) {
    for (int i = 0; i < nlayers; ++i) {
      int i1 = i + NGHOST;
      *(buf++) = tau_(b, i1);
      *(buf++) = ssa_(b, i1);
      for (int j = 0; j <= npmom; ++j) {
        *(buf++) = pmom_(b, i1, j);
      }
    }
  }
}

void RadiationBand::unpackSpectralProperties(int b, void *arg) {
  int nblocks = pexv->GetGroupSize();
  int nlayers = GetNumLayers();
  int npmom = GetNumPhaseMoments();

  auto buf = pexv->recv_buffer[1].data() + b * nlayers * (npmom + 3);

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
}
