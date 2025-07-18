// C/C++ headers
#include <iomanip>
#include <sstream>

// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/hydro/hydro.hpp>

// snap
#include <snap/eos/ideal_moist.hpp>

// canoe
#include <checks.hpp>
#include <impl.hpp>
#include <interface/eos.hpp>

// snap
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

void Decomposition::ChangeToAnomaly(AthenaArray<Real> &w, int kl, int ku,
                                    int jl, int ju) {
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;

  // positive in the x-increasing direction
  Real grav = pmb->phydro->hsrc.GetG1();
  Real gammad = pmb->peos->GetGamma();
  Real Rd = get_rd();

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

  // density
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = is + 1; i <= ie; ++i) {
        dsf_(k, j, i) =
            (w(IPR, k, j, i) - w(IPR, k, j, i - 1)) / (pco->dx1f(i) * grav);
      }

  auto pbot = ExchangeUtils::find_bot_neighbor(pmb);
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
  auto peos = pmb->pimpl->peos;
  auto cv_ratio_m1 = peos->cv_ratio_m1.accessor<Real, 1>();
  auto inv_mu_ratio_m1 = peos->inv_mu_ratio_m1.accessor<Real, 1>();
  int nvapor = peos->pthermo->options.vapor_ids().size() - 1;

  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j) {
      // adiabatic to bottom boundary
      // calculate local polytropic index
      Real fsig = 1., feps = 1.;
      for (int n = 1; n <= nvapor; ++n) {
        fsig += w(n, k, j, is) * cv_ratio_m1[n - 1];
        feps += w(n, k, j, is) * inv_mu_ratio_m1[n - 1];
      }
      Real gammas = 1. + (gammad - 1.) * feps / fsig;
      dsf_(k, j, is) = w(IDN, k, j, is) *
                       pow(psf_(k, j, is) / w(IPR, k, j, is), 1. / gammas);

      // isothermal to top boundary
      Real RdTv = w(IPR, k, j, ie) / w(IDN, k, j, ie);
      dsf_(k, j, ie + 1) =
          w(IDN, k, j, ie) * exp(grav * pco->dx1f(ie) / (2 * RdTv));

      // 1. change density and pressure (including ghost cells)
      for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
        // interpolate hydrostatic pressure, prevent divided by zero
        if (fabs(psf_(k, j, i) - psf_(k, j, i + 1)) < 1.E-6)
          psv_(k, j, i) = 0.5 * (psf_(k, j, i) + psf_(k, j, i + 1));
        else
          psv_(k, j, i) = (psf_(k, j, i) - psf_(k, j, i + 1)) /
                          log(psf_(k, j, i) / psf_(k, j, i + 1));

        // change pressure to pertubation
        w(IPR, k, j, i) -= psv_(k, j, i);
      }

      for (int i = is + 1; i <= ie - 1; ++i) {
        if (fabs(dsf_(k, j, i) - dsf_(k, j, i + 1)) < 1.E-6)
          dsv_(k, j, i) = 0.5 * (dsf_(k, j, i) + dsf_(k, j, i + 1));
        else
          dsv_(k, j, i) = (dsf_(k, j, i) - dsf_(k, j, i + 1)) /
                          log(dsf_(k, j, i) / dsf_(k, j, i + 1));

        // change density to anomaly
        // w(IDN, k, j, i) -= dsv_(k, j, i);
      }

      /* 2. adjust bottom boundary condition
      for (int i = is - NGHOST; i <= is; ++i) {
        dsv_(k, j, i) = w(IDN, k, j, i);
        w(IDN, k, j, i) = 0.;
      }

      if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::outflow) {
        for (int i = is - NGHOST; i < is; ++i) {
          w(IPR, k, j, i) = w(IPR, k, j, is);
        }
      }

      // 3. adjust top boundary condition
      for (int i = ie; i <= ie + NGHOST; ++i) {
        dsv_(k, j, i) = w(IDN, k, j, i);
        w(IDN, k, j, i) = 0.;
      }

      if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::outflow) {
        for (int i = ie + 1; i <= ie + NGHOST; ++i) {
          w(IPR, k, j, i) = w(IPR, k, j, ie);
        }
      }*/
    }

  /* debug
  if (Globals::my_rank == 0) {
    //int km = (kl + ku)/2;
    //int jm = (jl + ju)/2;
    int km = kl;
    int jm = jl;
    std::cout << "my.gid = " << pmb->gid << std::endl;
    std::cout << "bblock.gid = " << bblock.snb.gid << std::endl;
    std::cout << "===== k = " << km << " j = " << jm << std::endl;
    for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
      if (i == is)
        std::cout << "-------- ";
      if (i == 0)
        std::cout << "i = " << "-1/2 ";
      else if (i == 1)
        std::cout << "i = " << "+1/2 ";
      else
        std::cout << "i = " << i-1 << "+1/2 ";
      std::cout << "psf = " << psf_(km,jm,i) << ", ";
      std::cout << "dsf = " << dsf_(km,jm,i) << std::endl;
      std::cout << "i = " << i  << "    ";
      std::cout << " pre = " << w(IPR,km,jm,i)
                << " den = " << w(IDN,km,jm,i) << std::endl;
      if (i == ie)
        std::cout << "-------- ";
      if (i == ie + NGHOST) {
        std::cout << "i = " << i+1 << "+1/2 ";
        std::cout << "psf = " << psf_(km,jm,i+1) << ", ";
        std::cout << "dsf = " << dsf_(km,jm,i+1) << std::endl;
      }
    }
    std::cout << "==========" << std::endl;
  }*/

  // finish send top pressure
#ifdef MPI_PARALLEL
  MPI_Status status;
  if (has_bot_neighbor && (bblock.snb.rank != Globals::my_rank))
    MPI_Wait(&req_send_bot_, &status);
#endif
}

void Decomposition::RestoreFromAnomaly(AthenaArray<Real> &w,
                                       AthenaArray<Real> &wl,
                                       AthenaArray<Real> &wr, int k, int j,
                                       int il, int iu) {
  MeshBlock *pmb = pmy_block_;
  Hydro *phydro = pmb->phydro;

  Real Rd = get_rd();
  Real grav = phydro->hsrc.GetG1();
  int is = pmb->is, ie = pmb->ie;
  if (grav == 0.) return;

  for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
    w(IPR, k, j, i) += psv_(k, j, i);
    // w(IDN, k, j, i) += dsv_(k, j, i);
  }

  for (int i = il; i <= iu; ++i) {
    wr(IPR, i) += psf_(k, j, i);
    if (wr(IPR, i) < 0.) wr(IPR, i) = psf_(k, j, i);

    // wr(IDN, i) += dsf_(k, j, i);
    // if (wr(IDN, i) < 0.) wr(IDN, i) = dsf_(k, j, i);

    wl(IPR, i) += psf_(k, j, i);
    if (wl(IPR, i) < 0.) wl(IPR, i) = psf_(k, j, i);

    // wl(IDN, i) += dsf_(k, j, i);
    // if (wl(IDN, i) < 0.) wl(IDN, i) = dsf_(k, j, i);
  }

  check_decomposition(wl, wr, k, j, il, iu);
}
