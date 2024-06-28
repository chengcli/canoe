// C/C++ headers
#include <iomanip>
#include <sstream>

// Athena++ headers
#include <athena/coordinates/coordinates.hpp>
#include <athena/eos/eos.hpp>
#include <athena/hydro/hydro.hpp>

// utils
#include <utils/ndarrays.hpp>

// canoe
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
        psf(k, j, i + 1) = psf(k, j, i) - grav * w(IDN, k, j, i) * pco->dx1f(i);
}

inline void IntegrateDownwards(AthenaArray<Real> &psf,
                               AthenaArray<Real> const &w, Coordinates *pco,
                               Real grav, int kl, int ku, int jl, int ju,
                               int il, int iu) {
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = iu; i >= il; --i)
        psf(k, j, i) = psf(k, j, i + 1) + grav * w(IDN, k, j, i) * pco->dx1f(i);
}

void Decomposition::ChangeToPerturbation(AthenaArray<Real> &w, int kl, int ku,
                                         int jl, int ju) {
  // Need to integrate upwards
  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;

  Real grav = -pmb->phydro->hsrc.GetG1();  // positive downward pointing
  Real Rd = Thermodynamics::GetInstance()->GetRd();
  int is = pmb->is, ie = pmb->ie;
  if (grav == 0.) return;

  std::stringstream msg;
  if (NGHOST < 3) {
    msg << "### FATAL ERROR in function [Hydro::DecomposePressure]" << std::endl
        << "number of ghost cells (NGHOST) must be at least 3" << std::endl;
    ATHENA_ERROR(msg);
  }

  FindNeighbors();

  if (!has_top_neighbor) {
    // isothermal extrapolation
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j) {
        Real RdTv = w(IPR, k, j, ie) / w(IDN, k, j, ie);
        psf_(k, j, ie + 1) =
            w(IPR, k, j, ie) * exp(-grav * pco->dx1f(ie) / (2 * RdTv));
      }
  } else {
    RecvBuffer(psf_, kl, ku, jl, ju, ie + 1, ie + NGHOST + 1, tblock);
  }
  IntegrateDownwards(psf_, w, pco, grav, kl, ku, jl, ju, is, ie);

  // populate ghost cells
  SendBuffer(psf_, kl, ku, jl, ju);

  // boundary condition
  if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::reflect) {
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = 1; i <= NGHOST; ++i)
          psf_(k, j, is - i) = psf_(k, j, is + i);
  } else if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::outflow) {
    IntegrateDownwards(psf_, w, pco, 0., kl, ku, jl, ju, is - NGHOST, is - 1);
  }

  if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::reflect) {
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = 1; i <= NGHOST; ++i)
          psf_(k, j, ie + i + 1) = psf_(k, j, ie + 1 - i);
  } else if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::outflow) {
    IntegrateUpwards(psf_, w, pco, 0., kl, ku, jl, ju, ie + 1, ie + NGHOST);
  }

  if (has_bot_neighbor)
    RecvBuffer(psf_, kl, ku, jl, ju, is - NGHOST, is, bblock);

  // decompose pressure and density
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j) {
      // 1. change density and pressure (including ghost cells)
      for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
        // save pressure and density
        pres_(k, j, i) = w(IPR, k, j, i);
        // dens_(k, j, i) = w(IDN, k, j, i);

        // interpolate hydrostatic pressure, prevent divided by zero
        Real psv;
        if (fabs(psf_(k, j, i) - psf_(k, j, i + 1)) < 1.E-6) {
          psv = (psf_(k, j, i) + psf_(k, j, i + 1)) / 2.;
        } else {
          psv = (psf_(k, j, i) - psf_(k, j, i + 1)) /
                log(psf_(k, j, i) / psf_(k, j, i + 1));
        }

        // change pressure/density to pertubation quantities
        w(IPR, k, j, i) -= psv;
      }
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
      std::cout << "psf = " << psf_(km,jm,i) << std::endl;
      std::cout << "i = " << i  << "    ";
      std::cout << " pre = " << w(IPR,km,jm,i)
                << " den = " << w(IDN,km,jm,i) << std::endl;
      if (i == ie)
        std::cout << "-------- ";
      if (i == ie + NGHOST) {
        std::cout << "i = " << i+1 << "+1/2 ";
        std::cout << "psf = " << psf_(km,jm,i+1) << std::endl;
      }
    }
    std::cout << "==========" << std::endl;
  }*/

  // finish send top pressure
  WaitToFinishSend();
}

void Decomposition::RestoreFromPerturbation(AthenaArray<Real> &w,
                                            AthenaArray<Real> &wl,
                                            AthenaArray<Real> &wr, int k, int j,
                                            int il, int iu) {
  MeshBlock *pmb = pmy_block_;
  Hydro *phydro = pmb->phydro;
  auto pthermo = Thermodynamics::GetInstance();
  Real Rd = pthermo->GetRd();
  Real gammad = pmb->peos->GetGamma();
  int is = pmb->is, ie = pmb->ie;
  if (phydro->hsrc.GetG1() == 0.) return;

  for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
    w(IPR, k, j, i) = pres_(k, j, i);
    // w(IDN, k, j, i) = dens_(k, j, i);
  }

  for (int i = il; i <= iu; ++i) {
    wr(IPR, i) += psf_(k, j, i);
    if (wr(IPR, i) < 0.) wr(IPR, i) = psf_(k, j, i);

    wl(IPR, i) += psf_(k, j, i);
    if (wl(IPR, i) < 0.) wl(IPR, i) = psf_(k, j, i);

    wr(IDN, i) =
        w(IDN, k, j, i) * pow(wr(IPR, i) / w(IPR, k, j, i), 1. / gammad);
    wl(IDN, i) = w(IDN, k, j, i - 1) *
                 pow(wl(IPR, i) / w(IPR, k, j, i - 1), 1. / gammad);
  }

  /* debug
  //int km = (pmb->ks + pmb->ke)/2, jm = (pmb->js + pmb->je)/2;
  int km = pmb->ks-1, jm = pmb->js-1;
  if (Globals::my_rank == 1 && km == k & jm == j) {
    std::cout << "my.gid = " << pmb->gid << std::endl;
    std::cout << "bblock.gid = " << bblock.snb.gid << std::endl;
    std::cout << "===== k = " << km << " j = " << jm << std::endl;
    for (int i = il; i <= iu; ++i) {
      if (i == is)
        std::cout << "-------- ";
      if (i == 0)
        std::cout << "i = " << "-1/2 ";
      else if (i == 1)
        std::cout << "i = " << "+1/2 ";
      else
        std::cout << "i = " << i-1 << "+1/2 ";
      std::cout << "psf = " << psf_(km,jm,i)
                << " wl(IPR) = " << wl(IPR,i)
                << " wr(IPR) = " << wr(IPR,i)
                << " wl(IDN) = " << wl(IDN,i)
                << " wr(IDN) = " << wr(IDN,i)
                << std::endl;
      std::cout << "i = " << i  << "    ";
      std::cout << " pre = " << w(IPR,km,jm,i) << " den = " << w(IDN,km,jm,i) <<
  std::endl; if (i == ie) std::cout << "-------- "; if (i == ie + NGHOST) {
        std::cout << "i = " << i+1 << "+1/2 ";
        std::cout << "psf = " << psf_(km,jm,i+1) << std::endl;
      }
    }
    std::cout << "==========" << std::endl;
  }*/
}
