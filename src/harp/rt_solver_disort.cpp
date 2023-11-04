//! \file rt_solver_disort.cpp
//! \brief Call DISORT to perform radiative transfer calculation

// C/C++
#include <cmath>
#include <iostream>

// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/mesh/mesh.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// canoe
#include <configure.hpp>
#include <impl.hpp>

// astro
#include <astro/celestrial_body.hpp>

// climath
#include <climath/special.h>

// harp
#include "radiation.hpp"
#include "rt_solvers.hpp"

#ifdef RT_DISORT

void RadiationBand::RTSolverDisort::CalBandFlux(MeshBlock const *pmb, int k,
                                                int j, int il, int iu) {
  auto pcoord = pmb->pcoord;
  auto planet = pmb->pimpl->planet;
  auto prad = pmb->pimpl->prad;

  auto &bflxup = pmy_band_->bflxup;
  auto &bflxdn = pmy_band_->bflxdn;
  auto &bpmom = pmy_band_->bpmom;

  auto &wmin = pmy_band_->wrange_.first;
  auto &wmax = pmy_band_->wrange_.second;
  auto &spec = pmy_band_->pgrid_->spec;

  auto &flxup = pmy_band_->flxup_;
  auto &flxdn = pmy_band_->flxdn_;

  Direction ray;
  Real dist_au;
  Real time = pmb->pmy_mesh->time;

  if (pmy_band_->TestFlag(RadiationFlags::Dynamic)) {
    ray = planet->ParentZenithAngle(time, pcoord->x2v(j), pcoord->x3v(k));
    dist_au = planet->ParentDistanceInAu(time);
  } else {
    ray = prad->GetRayInput(0);
    dist_au = pmb->pimpl->GetDistanceInAu();
  }

  if (ds_.flag.ibcnd != 0) {
    throw ValueError("RTSolverDisort::CalRadtranFlux", "ibcnd", ds_.flag.ibcnd,
                     0);
  }
  pmy_band_->Regroup(pmb, X1DIR);

  // transfer temperature
  if (ds_.flag.planck) {
    pmy_band_->PackTemperature();
    pmy_band_->Transfer(pmb, 0);
    pmy_band_->UnpackTemperature(&ds_);
  }

  // for (int i = 0; i <= ds_.nlyr; ++i)
  //   std::cout << ds_.temper[i] << std::endl;

  ds_.bc.umu0 = ray.mu > 1.E-3 ? ray.mu : 1.E-3;
  ds_.bc.phi0 = ray.phi;
  if (ds_.flag.planck) {
    ds_.bc.btemp = ds_.temper[ds_.nlyr];
    ds_.bc.ttemp = ds_.temper[0];
  }

  // reset flx of this column
  for (int i = il; i <= iu; ++i) bflxup(k, j, i) = bflxdn(k, j, i) = 0.;

  AthenaArray<Real> farea(iu + 1), vol(iu + 1);
  pcoord->Face1Area(k, j, il, iu, farea);
  pcoord->CellVolume(k, j, il, iu, vol);

  if (pmy_band_->TestFlag(RadiationFlags::CorrelatedK)) {
    // stellar source function
    if (pmy_band_->TestFlag(RadiationFlags::Star))
      ds_.bc.fbeam = planet->ParentInsolationFlux(wmin, wmax, dist_au);
    // planck source function
    ds_.wvnmlo = wmin;
    ds_.wvnmhi = wmax;
  }

  int rank = pmy_band_->GetRankInGroup();

  // loop over spectral grids in the band
  for (int n = 0; n < pmy_band_->GetNumSpecGrids(); ++n) {
    if (!(pmy_band_->TestFlag(RadiationFlags::CorrelatedK))) {
      // stellar source function
      if (pmy_band_->TestFlag(RadiationFlags::Star))
        ds_.bc.fbeam =
            planet->ParentInsolationFlux(spec[n].wav1, spec[n].wav2, dist_au);
      // planck source function
      ds_.wvnmlo = spec[n].wav1;
      ds_.wvnmhi = spec[n].wav2;
    }

    // transfer spectral grid data
    pmy_band_->PackSpectralGrid(n);
    pmy_band_->Transfer(pmb, 1);
    pmy_band_->UnpackSpectralGrid(&ds_);

    // run disort
    c_disort(&ds_, &ds_out_);

    /// \note Counting Index
    ///
    /// Example, il = 0, iu = 2, ds_.nlyr = 6, partition in to 3 blocks
    /// face id   -> 0 - 1 - 2 - 3 - 4 - 5 - 6
    /// cell id   -> | 0 | 1 | 2 | 3 | 4 | 5 |
    /// disort id -> 6 - 5 - 4 - 3 - 2 - 1 - 0
    /// blocks    -> ---------       *       *
    ///           ->  r = 0  *       *       *
    ///           ->         ---------       *
    ///           ->           r = 1 *       *
    ///           ->                 ---------
    ///           ->                   r = 2
    /// block r = 0 gets, 6 - 5 - 4
    /// block r = 1 gets, 4 - 3 - 2
    /// block r = 2 gets, 2 - 1 - 0
    /// accumulate flux from lines
    for (int i = il; i <= iu; ++i) {
      int m = ds_.nlyr - (rank * (iu - il) + i - il);
      //! \bug does not work for spherical geometry, need to scale area using
      //! farea(il)/farea(i)
      // flux up
      flxup(n, i) = ds_out_.rad[m].flup;

      //! \bug does not work for spherical geomtry, need to scale area using
      //! farea(il)/farea(i)
      // flux down
      flxdn(n, i) = ds_out_.rad[m].rfldir + ds_out_.rad[m].rfldn;

      bflxup(k, j, i) += spec[n].wght * flxup(n, i);
      bflxdn(k, j, i) += spec[n].wght * flxdn(n, i);
    }

    // spherical correction by XIZ
    // xiz 2022 flux scaling so that the heating rate is the same as the
    // plane-parallel scheme volheating scaling: first calculate flux divergence
    // from DISORT using Plane-parallel in a cell then mulitpled by the cell
    // volume divided by dx1f then solve for F using F1*S1-F2*S2 = volheating
    // the top fluxes are the still the same as the plane-paralell values
    Real volh, bflxup1 = bflxup(k, j, iu), bflxdn1 = bflxdn(k, j, iu);
    for (int i = iu - 1; i >= il; --i) {
      // upward
      volh = (bflxup1 - bflxup(k, j, i)) / pcoord->dx1f(i) * vol(i);
      bflxup1 = bflxup(k, j, i);
      bflxup(k, j, i) = (bflxup(k, j, i + 1) * farea(i + 1) - volh) / farea(i);

      // downward
      volh = (bflxdn1 - bflxdn(k, j, i)) / pcoord->dx1f(i) * vol(i);
      bflxdn1 = bflxdn(k, j, i);
      bflxdn(k, j, i) = (bflxdn(k, j, i + 1) * farea(i + 1) - volh) / farea(i);
    }
  }
}

void RadiationBand::RTSolverDisort::CalBandRadiance(MeshBlock const *pmb, int k,
                                                    int j, int il, int iu) {
  throw NotImplementedError("RTSolverDisort::CalBandRadiance");
}

#endif  // RT_DISORT
