// C/C++ headers
#include <iomanip>
#include <sstream>

// Athena++ headers
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/hydro/hydro.hpp>

// utils
#include <utils/ndarrays.hpp>

// exchanger
#include <exchanger/exchanger.hpp>

// canoe
#include <checks.hpp>
#include <impl.hpp>

// snap
#include "../thermodynamics/thermodynamics.hpp"
#include "decomposition.hpp"

inline void IntegrateUpwards(AthenaArray<Real> &psf, AthenaArray<Real> const &w,
                             Coordinates *pco, Real grav, int kl, int ku,
                             int jl, int ju, int il, int iu) {
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = il; i <= iu; ++i)
        psf(k, j, i + 1) = psf(k, j, i) + grav * w(IDN, k, j, i) * pco->dx1f(i);
}

inline void IntegrateDownwards(AthenaArray<Real> &psf,
                               AthenaArray<Real> const &w, Coordinates *pco,
                               Real grav, int kl, int ku, int jl, int ju,
                               int il, int iu) {
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = iu; i >= il; --i)
        psf(k, j, i) = psf(k, j, i + 1) - grav * w(IDN, k, j, i) * pco->dx1f(i);
}

void Decomposition::ChangeToTemperature(AthenaArray<Real> &w, int kl, int ku,
                                        int jl, int ju) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;
  auto pthermo = Thermodynamics::GetInstance();

  // positive in the x-increasing direction
  Real grav = pmb->phydro->hsrc.GetG1();
  Real Rd = pthermo->GetRd();

  int is = pmb->is, ie = pmb->ie;
  if (grav == 0.) return;

  FindNeighbors();

  if (has_top_neighbor) {
    RecvFromTop(psf_, kl, ku, jl, ju);
  } else {
    // isothermal extrapolation to find the pressure at top boundary
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j) {
        Real RdTv = w(IPR, k, j, ie) / w(IDN, k, j, ie);
        psf_(k, j, ie + 1) =
            w(IPR, k, j, ie) * exp(grav * pco->dx1f(ie) / (2 * RdTv));
      }
  }
  IntegrateDownwards(psf_, w, pco, grav, kl, ku, jl, ju, is, ie);

  // populate ghost cells
  if (has_bot_neighbor) SendToBottom(psf_, kl, ku, jl, ju);

  // boundary condition
  if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::reflect) {
    IntegrateDownwards(psf_, w, pco, -grav, kl, ku, jl, ju, is - NGHOST,
                       is - 1);
  } else {  // block boundary
    IntegrateDownwards(psf_, w, pco, grav, kl, ku, jl, ju, is - NGHOST, is - 1);
  }

  if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::reflect) {
    IntegrateUpwards(psf_, w, pco, -grav, kl, ku, jl, ju, ie + 1, ie + NGHOST);
  } else {  // block boundary
    IntegrateUpwards(psf_, w, pco, grav, kl, ku, jl, ju, ie + 1, ie + NGHOST);
  }

  // decompose pressure and density
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j) {
      // 1. change density and pressure (including ghost cells)
      for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
        // interpolate hydrostatic pressure, prevent divided by zero
        if (fabs(psf_(k, j, i) - psf_(k, j, i + 1)) < 1.E-6)
          psv_(k, j, i) = 0.5 * (psf_(k, j, i) + psf_(k, j, i + 1));
        else
          psv_(k, j, i) = (psf_(k, j, i) - psf_(k, j, i + 1)) /
                          log(psf_(k, j, i) / psf_(k, j, i + 1));

        // calculate RovRd
        Real feps = 1.;
        for (int n = 1; n <= NVAPOR; ++n)
          feps += w(n, k, j, i) * (pthermo->GetInvMuRatio(n) - 1.);
        gamma_(k, j, i) = Rd * feps;

        // change density to temperature
        w(IDN, k, j, i) = w(IPR, k, j, i) / (w(IDN, k, j, i) * gamma_(k, j, i));

        // change pressure to pertubation
        w(IPR, k, j, i) -= psv_(k, j, i);
      }

      // 2. adjust bottom boundary condition
      if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::outflow) {
        for (int i = is - NGHOST; i < is; ++i) {
          w(IDN, k, j, i) = w(IDN, k, j, is);
          w(IPR, k, j, i) = w(IPR, k, j, is);
        }
      }

      // 3. adjust top boundary condition
      if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::outflow) {
        for (int i = ie + 1; i <= ie + NGHOST; ++i) {
          w(IDN, k, j, i) = w(IDN, k, j, ie);
          w(IPR, k, j, i) = w(IPR, k, j, ie);
        }
      }
    }

    // finish send top pressure
#ifdef MPI_PARALLEL
  MPI_Status status;
  if (has_bot_neighbor && (bblock.snb.rank != Globals::my_rank))
    MPI_Wait(&req_send_bot_, &status);
#endif
}

void Decomposition::RestoreFromTemperature(AthenaArray<Real> &w,
                                           AthenaArray<Real> &wl,
                                           AthenaArray<Real> &wr, int k, int j,
                                           int il, int iu) {
  MeshBlock *pmb = pmy_block_;
  Hydro *phydro = pmb->phydro;
  auto pthermo = Thermodynamics::GetInstance();

  Real Rd = pthermo->GetRd();
  Real grav = phydro->hsrc.GetG1();
  int is = pmb->is, ie = pmb->ie;
  if (grav == 0.) return;

  for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
    w(IPR, k, j, i) += psv_(k, j, i);
    w(IDN, k, j, i) = w(IPR, k, j, i) / (w(IDN, k, j, i) * gamma_(k, j, i));
  }

  for (int i = il; i <= iu; ++i) {
    wr(IPR, i) += psf_(k, j, i);
    if (wr(IPR, i) < 0.) wr(IPR, i) = psf_(k, j, i);
    Real gamma = sqrt(gamma_(k, j, i) * gamma_(k, j, i - 1));
    wr(IDN, i) = wr(IPR, i) / (wr(IDN, i) * gamma);

    wl(IPR, i) += psf_(k, j, i);
    if (wl(IPR, i) < 0.) wl(IPR, i) = psf_(k, j, i);
    wl(IDN, i) = wl(IPR, i) / (wl(IDN, i) * gamma);
  }

  check_decomposition(wl, wr, k, j, il, iu);
}
