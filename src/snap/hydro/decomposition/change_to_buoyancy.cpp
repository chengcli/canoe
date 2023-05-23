// C/C++ headers
#include <iomanip>
#include <sstream>

// Athena++ headers
#include <coordinates/coordinates.hpp>

// cliuitls
#include <cliutils/ndarrays.hpp>
#include <cliutils/stride_iterator.hpp>

// canoe headers
#include "../../mesh/meshblock_impl.hpp"
#include "../../thermodynamics/thermodynamics.hpp"
#include "decomposition.hpp"

inline void IntegrateDownwards(AthenaArray<Real> &psf,
                               AthenaArray<Real> const &w, Coordinates *pco,
                               Real grav, int kl, int ku, int jl, int ju,
                               int il, int iu) {
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = iu; i >= il; --i)
        psf(k, j, i) = psf(k, j, i + 1) + grav * w(IDN, k, j, i) * pco->dx1f(i);
}

void Decomposition::ChangeToBuoyancy(AthenaArray<Real> &w, int kl, int ku,
                                     int jl, int ju) {
  MeshBlock *pmb = pmy_hydro->pmy_block;
  Coordinates *pco = pmb->pcoord;
  Thermodynamics *pthermo = pmb->pimpl->pthermo;

  Real grav = -pmy_hydro->hsrc.GetG1();  // positive downward pointing

  int is = pmb->is, ie = pmb->ie;

  if (grav == 0.) return;

  std::stringstream msg;
  if (NGHOST < 3) {
    msg << "### FATAL ERROR in function [Hydro::DecomposePressure]" << std::endl
        << "number of ghost cells (NGHOST) must be at least 3" << std::endl;
    ATHENA_ERROR(msg);
  }

  Real **w1;
  NewCArray(w1, 2, NHYDRO);
  FindNeighbors();

  if (has_top_neighbor) {
    RecvFromTop(psf_, kl, ku, jl, ju);
  } else {
    // adiabatic extrapolation
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j) {
        Real P1 = w(IPR, k, j, ie);
        Real T1 = pthermo->getTemp(w.at(k, j, ie));
        Real dz = pco->dx1f(ie);
        for (int n = 0; n < NHYDRO; ++n) w1[0][n] = w(n, k, j, ie);

        // adiabatic extrapolation for half a grid
        pthermo->ConstructAtmosphere(w1, T1, P1, grav, dz / 2., 2,
                                     Adiabat::reversible, 0.);
        psf_(k, j, ie + 1) = w1[1][IPR];
      }
  }
  IntegrateDownwards(psf_, w, pco, grav, kl, ku, jl, ju, is, ie);

  if (has_bot_neighbor) SendToBottom(psf_, kl, ku, jl, ju);

  // decompose pressure
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j) {
      // save pressure and density
      for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
        pres_(k, j, i) = w(IPR, k, j, i);
        dens_(k, j, i) = w(IDN, k, j, i);
      }
      // decompose pressure and density
      for (int i = is; i <= ie; ++i) {
        Real psv;
        // interpolate hydrostatic pressure, prevent divided by zero
        if (fabs(psf_(k, j, i) - psf_(k, j, i + 1)) < 1.E-6)
          psv = (psf_(k, j, i) + psf_(k, j, i + 1)) / 2.;
        else
          psv = (psf_(k, j, i) - psf_(k, j, i + 1)) /
                log(psf_(k, j, i) / psf_(k, j, i + 1));
        w(IPR, k, j, i) -= psv;

        w(IDN, k, j, i) = (pres_(k, j, i - 1) - pres_(k, j, i + 1)) /
                              (pco->x1v(i + 1) - pco->x1v(i - 1)) -
                          dens_(k, j, i) * grav;
      }

      // fix boundary condition
      if (pmb->pbval->block_bcs[inner_x1] != BoundaryFlag::block)
        w(IDN, k, j, is) = (pres_(k, j, is) - pres_(k, j, is + 1)) /
                               (pco->x1v(is + 1) - pco->x1v(is)) -
                           sqrt(dens_(k, j, is) * dens_(k, j, is + 1)) * grav;

      if (pmb->pbval->block_bcs[outer_x1] != BoundaryFlag::block)
        w(IDN, k, j, ie) = (pres_(k, j, ie - 1) - pres_(k, j, ie)) /
                               (pco->x1v(ie) - pco->x1v(ie - 1)) -
                           sqrt(dens_(k, j, ie) * dens_(k, j, ie - 1)) * grav;
    }

  SyncNewVariables(w, kl, ku, jl, ju);
  ApplyHydroBoundary(w, psf_, kl, ku, jl, ju);

  /* debug
  if (Globals::my_rank == 0) {
    //int km = (kl + ku)/2;
    //int jm = (jl + ju)/2;
    int km = ku;
    int jm = ju;
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
      std::cout << "psf = " << std::setprecision(8) << psf_(km,jm,i) <<
  std::endl; std::cout << "i = " << i  << "    "; std::cout << " pre = " <<
  w(IPR,km,jm,i) << " den = " << w(IDN,km,jm,i) << std::endl; if (i == ie)
        std::cout << "-------- ";
      if (i == ie + NGHOST) {
        std::cout << "i = " << i+1 << "+1/2 ";
        std::cout << "psf = " << psf_(km,jm,i+1) << std::endl;
      }
    }
    std::cout << "==========" << std::endl;
  }*/

  FreeCArray(w1);
  // finish MPI communication
  WaitToFinishSend();
  WaitToFinishSync(w, kl, ku, jl, ju);
}

void Decomposition::RestoreFromBuoyancy(AthenaArray<Real> &w,
                                        AthenaArray<Real> &wl,
                                        AthenaArray<Real> &wr, int k, int j,
                                        int il, int iu) {
  MeshBlock *pmb = pmy_hydro->pmy_block;
  Coordinates *pco = pmb->pcoord;
  Thermodynamics *pthermo = pmb->pimpl->pthermo;
  Real grav = -pmy_hydro->hsrc.GetG1();  // positive downward pointing
  int is = pmb->is, ie = pmb->ie;
  if (grav == 0.) return;

  for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
    w(IPR, k, j, i) = pres_(k, j, i);
    w(IDN, k, j, i) = dens_(k, j, i);
  }

  Real mdpdz;
  for (int i = il; i <= iu; ++i) {
    wr(IPR, i) += psf_(k, j, i);
    if (wr(IPR, i) < 0.) wr(IPR, i) = psf_(k, j, i);

    wl(IPR, i + 1) += psf_(k, j, i + 1);
    if (wl(IPR, i + 1) < 0.) wl(IPR, i + 1) = psf_(k, j, i + 1);

    mdpdz = (w(IPR, k, j, i - 1) - w(IPR, k, j, i)) /
            (pco->x1v(i) - pco->x1v(i - 1));
    wr(IDN, i) = (mdpdz - wr(IDN, i)) / grav;
    if (wr(IDN, i) < 0.) wr(IDN, i) = mdpdz / grav;

    wl(IDN, i + 1) = (mdpdz - wl(IDN, i + 1)) / grav;
    if (wl(IDN, i + 1) < 0.) wl(IDN, i + 1) = mdpdz / grav;
  }

  // fix boundary condition
  Real **w1;
  NewCArray(w1, 2, NHYDRO);

  // adiabatic extrapolation for a grid
  Real P1 = w(IPR, k, j, is);
  Real T1 = pthermo->getTemp(w.at(k, j, is));
  Real dz = pco->dx1f(is);
  for (int n = 0; n < NHYDRO; ++n) w1[0][n] = w(n, k, j, is);
  pthermo->ConstructAtmosphere(w1, T1, P1, grav, -dz, 2, Adiabat::reversible,
                               0.);

  mdpdz = (w1[1][IPR] - w1[0][IPR]) / dz;
  if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::reflect) {
    for (int i = il + 1; i < is; ++i) {
      wl(IDN, i) = wl(IDN, 2 * is - i);
      wr(IDN, i) = wr(IDN, 2 * is - i);
    }
    wl(IDN, is) = (mdpdz - wl(IDN, is)) / grav;
    wr(IDN, is) = (mdpdz - wr(IDN, is)) / grav;
  } else if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::outflow) {
    for (int i = il + 1; i < is; ++i)
      wl(IDN, i) = wr(IDN, i) = sqrt(w(IDN, k, j, i) * w(IDN, k, j, i - 1));
    wl(IDN, is) = (mdpdz - wl(IDN, is)) / grav;
    wr(IDN, is) = (mdpdz - wr(IDN, is)) / grav;
  }

  // adiabatic extrapolation for a grid
  P1 = w(IPR, k, j, ie);
  T1 = pthermo->getTemp(w.at(k, j, ie));
  dz = pco->dx1f(ie);
  for (int n = 0; n < NHYDRO; ++n) w1[0][n] = w(n, k, j, ie);
  pthermo->ConstructAtmosphere(w1, T1, P1, grav, dz, 2, Adiabat::reversible,
                               0.);

  mdpdz = (w1[0][IPR] - w1[1][IPR]) / dz;
  if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::reflect) {
    for (int i = ie + 2; i <= iu + 1; ++i) {
      wl(IDN, i) = wl(IDN, 2 * ie - i + 2);
      wr(IDN, i) = wr(IDN, 2 * ie - i + 2);
    }
    wl(IDN, ie + 1) = (mdpdz - wl(IDN, ie + 1)) / grav;
    wr(IDN, ie + 1) = (mdpdz - wr(IDN, ie + 1)) / grav;
  } else if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::outflow) {
    for (int i = ie + 2; i <= iu + 1; ++i)
      wl(IDN, i) = wr(IDN, i) = sqrt(w(IDN, k, j, i) * w(IDN, k, j, i - 1));
    wl(IDN, ie + 1) = (mdpdz - wl(IDN, ie + 1)) / grav;
    wr(IDN, ie + 1) = (mdpdz - wr(IDN, ie + 1)) / grav;
  }

  FreeCArray(w1);

  /* debug
  int km = (pmb->ks + pmb->ke)/2, jm = (pmb->js + pmb->je)/2;
  if (Globals::my_rank == 0 && km == k & jm == j) {
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
