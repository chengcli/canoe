//! \file rt_solver_disort.cpp
//! \brief Call DISORT to perform radiative transfer calculation

// C/C++
#include <cmath>
#include <iostream>

// external
#include <yaml-cpp/yaml.h>

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

// harp
#include "radiation.hpp"
#include "rt_solvers.hpp"

#ifdef RT_DISORT

RadiationBand::RTSolverDisort::RTSolverDisort(RadiationBand *pmy_band,
                                              YAML::Node const &rad)
    : RTSolver(pmy_band, "Disort") {
  if (rad["Disort"]) setFlagsFromNode(rad["Disort"]);
}

//! \todo update based on band outdir
void RadiationBand::RTSolverDisort::Resize(int nlyr, int nstr, int nuphi,
                                           int numu) {
  SetAtmosphereDimension(nlyr, nstr, nstr, nstr);
  SetIntensityDimension(nuphi, 1, numu);
  Finalize();

  Real utau = 0.;
  Real uphi = 0.;
  Real umu = 1.;

  SetUserAzimuthalAngle(&uphi, 1);
  SetUserOpticalDepth(&utau, 1);
  SetUserCosinePolarAngle(&umu, 1);
}

//! \note Counting Disort Index
//! Example, il = 0, iu = 2, ds_.nlyr = 6, partition in to 3 blocks
//! face id   -> 0 - 1 - 2 - 3 - 4 - 5 - 6
//! cell id   -> | 0 | 1 | 2 | 3 | 4 | 5 |
//! disort id -> 6 - 5 - 4 - 3 - 2 - 1 - 0
//! blocks    -> ---------       *       *
//!           ->  r = 0  *       *       *
//!           ->         ---------       *
//!           ->           r = 1 *       *
//!           ->                 ---------
//!           ->                   r = 2
//! block r = 0 gets, 6 - 5 - 4
//! block r = 1 gets, 4 - 3 - 2
//! block r = 2 gets, 2 - 1 - 0

void RadiationBand::RTSolverDisort::Prepare(MeshBlock const *pmb, int k,
                                            int j) {
  auto &wmin = pmy_band_->wrange_.first;
  auto &wmax = pmy_band_->wrange_.second;
  auto planet = pmb->pimpl->planet;
  auto pcoord = pmb->pcoord;

  Direction ray;
  Real dist_au;
  Real time = pmb->pmy_mesh->time;

  if (pmy_band_->TestFlag(RadiationFlags::Dynamic)) {
    ray = planet->ParentZenithAngle(time, pcoord->x2v(j), pcoord->x3v(k));
    dist_au = planet->ParentDistanceInAu(time);
  } else {
    ray = pmb->pimpl->prad->GetRayInput(0);
    dist_au = pmb->pimpl->GetDistanceInAu();
  }

  if (ds_.flag.ibcnd != 0) {
    throw ValueError("RTSolverDisort::CalRadtranFlux", "ibcnd", ds_.flag.ibcnd,
                     0);
  }

  pmy_band_->Regroup(pmb, X1DIR);
  myrank_in_column_ = pmy_band_->GetRankInGroup();

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

  if (pmy_band_->TestFlag(RadiationFlags::CorrelatedK)) {
    // stellar source function
    if (pmy_band_->TestFlag(RadiationFlags::Star))
      ds_.bc.fbeam = planet->ParentInsolationFlux(wmin, wmax, dist_au);
    // planck source function
    ds_.wvnmlo = wmin;
    ds_.wvnmhi = wmax;
  }

  pmb->pcoord->Face1Area(k, j, pmb->is, pmb->ie + 1, farea_);
  pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol_);

  Finalize();
}

void RadiationBand::RTSolverDisort::CalBandFlux(MeshBlock const *pmb, int k,
                                                int j, int il, int iu) {
  Real dist_au;
  if (pmy_band_->TestFlag(RadiationFlags::Dynamic)) {
    dist_au = pmb->pimpl->planet->ParentDistanceInAu(time);
  } else {
    dist_au = pmb->pimpl->GetDistanceInAu();
  }

  // loop over spectral grids in the band
  int b = 0;
  for (auto &spec : pmy_band_->pgrid_->spec) {
    if (!(pmy_band_->TestFlag(RadiationFlags::CorrelatedK))) {
      // stellar source function
      if (pmy_band_->TestFlag(RadiationFlags::Star)) {
        ds_.bc.fbeam = pmb->pimpl->planet->ParentInsolationFlux(
            spec.wav1, spec.wav2, dist_au);
      }
      // planck source function
      ds_.wvnmlo = spec.wav1;
      ds_.wvnmhi = spec.wav2;
    }

    // transfer spectral grid data
    pmy_band_->PackSpectralGrid(b);
    pmy_band_->Transfer(pmb, 1);
    pmy_band_->UnpackSpectralGrid(&ds_);

    // run disort
    c_disort(&ds_, &ds_out_);

    // add spectral bin flux
    addDisortFlux(pmb->pcoord, b++, k, j, il, iu);
  }
}

void RadiationBand::RTSolverDisort::addDisortFlux(Coordinates const *pcoord,
                                                  int n, int k, int j, int il,
                                                  int iu) {
  auto &bflxup = pmy_band_->bflxup;
  auto &bflxdn = pmy_band_->bflxdn;

  auto &flxup = pmy_band_->flxup_;
  auto &flxdn = pmy_band_->flxdn_;
  auto &spec = pmy_band_->pgrid_->spec;

  /// accumulate flux from spectral bins
  for (int i = il; i <= iu; ++i) {
    int m = ds_.nlyr - (myrank_in_column_ * (iu - il) + i - il);
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

  //! \note Spherical correction by XIZ
  //! xiz 2022 flux scaling so that the heating rate is the same as the
  //! plane-parallel scheme volheating scaling: first calculate flux divergence
  //! from DISORT using Plane-parallel in a cell then mulitpled by the cell
  //! volume divided by dx1f then solve for F using F1*S1-F2*S2 = volheating
  //! the top fluxes are the still the same as the plane-paralell values
  Real volh;
  Real bflxup1 = bflxup(k, j, iu);
  Real bflxdn1 = bflxdn(k, j, iu);

  for (int i = iu - 1; i >= il; --i) {
    // upward
    volh = (bflxup1 - bflxup(k, j, i)) / pcoord->dx1f(i) * vol_(i);
    bflxup1 = bflxup(k, j, i);
    bflxup(k, j, i) = (bflxup(k, j, i + 1) * farea_(i + 1) - volh) / farea_(i);

    // downward
    volh = (bflxdn1 - bflxdn(k, j, i)) / pcoord->dx1f(i) * vol_(i);
    bflxdn1 = bflxdn(k, j, i);
    bflxdn(k, j, i) = (bflxdn(k, j, i + 1) * farea_(i + 1) - volh) / farea_(i);
  }
}

void RadiationBand::RTSolverDisort::CalBandRadiance(MeshBlock const *pmb, int k,
                                                    int j) {
  if (ds_.flag.onlyfl) {
    throw RuntimeError("RTSolverDisort::CalBandRadiance",
                       "Only flux calculation is requested");
  }

  if (ds_.ntau != 1) {
    throw RuntimeError("RTSolverDisort::CalBandRadiance",
                       "Only toa radiance (ds.ntau = 1) is supported");
  }

  int nrays = ds_.nphi * ds_.ntau * ds_.numu;

  if (nrays != pmy_band_->GetNumOutgoingRays()) {
    throw RuntimeError("RTSolverDisort::CalBandRadiance",
                       "Number of outgoing rays does not match");
  }

  Real dist_au;
  if (pmy_band_->TestFlag(RadiationFlags::Dynamic)) {
    dist_au = pmb->pimpl->planet->ParentDistanceInAu(time);
  } else {
    dist_au = pmb->pimpl->GetDistanceInAu();
  }

  // loop over spectral grids in the band
  int b = 0;
  for (auto &spec : pmy_band_->pgrid_->spec) {
    if (!(pmy_band_->TestFlag(RadiationFlags::CorrelatedK))) {
      // stellar source function
      if (pmy_band_->TestFlag(RadiationFlags::Star)) {
        ds_.bc.fbeam = pmb->pimpl->planet->ParentInsolationFlux(
            spec.wav1, spec.wav2, dist_au);
      }
      // planck source function
      ds_.wvnmlo = spec.wav1;
      ds_.wvnmhi = spec.wav2;
    }

    // transfer spectral grid data
    pmy_band_->PackSpectralGrid(b);
    pmy_band_->Transfer(pmb, 1);
    pmy_band_->UnpackSpectralGrid(&ds_);

    // run disort
    c_disort(&ds_, &ds_out_);

    // add spectral bin radiance
    addDisortRadiance(pmb->pcoord, b++, k, j);
  }
}

void RadiationBand::RTSolverDisort::addDisortRadiance(Coordinates const *pcoord,
                                                      int b, int k, int j) {
  auto &btoa = pmy_band_->btoa;
  auto &spec = pmy_band_->pgrid_->spec;

  for (int c = 0; c < ds_.nphi * ds_.ntau * ds_.numu; ++c) {
    btoa(c, k, j) += spec[b].wght * ds_out_.uu[c];
  }
}

#endif  // RT_DISORT
