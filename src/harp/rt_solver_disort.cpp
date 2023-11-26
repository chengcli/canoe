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

// climath
#include <climath/interpolation.h>

// canoe
#include <configure.hpp>
#include <constants.hpp>
#include <impl.hpp>

// astro
#include <astro/celestrial_body.hpp>

// harp
#include "radiation.hpp"
#include "rt_solvers.hpp"

#ifdef RT_DISORT

std::map<std::string, bool> to_map_bool(YAML::Node const &node) {
  std::map<std::string, bool> flags;

  for (auto it = node.begin(); it != node.end(); ++it) {
    flags[it->first.as<std::string>()] = it->second.as<bool>();
  }

  return flags;
}

RadiationBand::RTSolverDisort::RTSolverDisort(RadiationBand *pmy_band,
                                              YAML::Node const &rad)
    : RTSolver(pmy_band, "Disort") {
  if (rad["Disort-flags"]) {
    SetFlags(to_map_bool(rad["Disort-flags"]));
  }

  // override disort planck flag
  if (pmy_band->TestFlag(RadiationFlags::ThermalEmission)) {
    ds_.flag.planck = true;
  }

  SetHeader("Disort solving band " + pmy_band_->GetName());
}

//! \todo update based on band outdir
void RadiationBand::RTSolverDisort::Resize(int nlyr, int nstr) {
  Unseal();

  auto &rayout = pmy_band_->rayOutput_;
  auto &&uphi_umu = RadiationHelper::get_direction_grids(rayout);

  SetAtmosphereDimension(nlyr, nstr, nstr, nstr);

  dir_dim_[0] = uphi_umu.second.size();  // umu
  dir_dim_[1] = uphi_umu.first.size();   // uphi
  dir_axis_.resize(dir_dim_[0] + dir_dim_[1]);

  SetIntensityDimension(dir_dim_[1], 1, dir_dim_[0]);
  Seal();
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

  Direction ray;
  Real dist_au = 1.;

  if (pmy_band_->TestFlag(RadiationFlags::StellarBeam)) {
    if (pmb == nullptr) {
      throw RuntimeError("RTSolverDisort::Prepare",
                         "MeshBlock must be allocated for stellar radiation");
    }
    Real time = pmb->pmy_mesh->time;
    auto planet = pmb->pimpl->planet;
    if (pmy_band_->TestFlag(RadiationFlags::TimeDependent)) {
      ray = planet->ParentZenithAngle(time, pmb->pcoord->x2v(j),
                                      pmb->pcoord->x3v(k));
      dist_au = planet->ParentDistanceInAu(time);
    } else if (pmb != nullptr) {
      ray = pmb->pimpl->prad->GetRayInput(0);
      dist_au = pmb->pimpl->GetDistanceInAu();
    }
  } else {
    if (pmy_band_->HasPar("umu0")) {
      ray.mu = pmy_band_->GetPar<Real>("umu0");
    }

    if (pmy_band_->HasPar("phi0")) {
      ray.phi = pmy_band_->GetPar<Real>("phi0");
    }

    if (pmy_band_->HasPar("dist_au")) {
      dist_au = pmy_band_->GetPar<Real>("dist_au");
    }
  }

  if (ds_.flag.ibcnd != 0) {
    throw ValueError("RTSolverDisort::CalRadtranFlux", "ibcnd", ds_.flag.ibcnd,
                     0);
  }

  if (pmb != nullptr) {
    pmy_band_->Regroup(pmb, X1DIR);
    myrank_in_column_ = pmy_band_->GetRankInGroup();
  } else {
    myrank_in_column_ = 0;
  }

  // transfer temperature
  if (ds_.flag.planck) {
    pmy_band_->PackTemperature();
    pmy_band_->Transfer(pmb, 0);
    pmy_band_->UnpackTemperature(&ds_);
  }

  ds_.bc.umu0 = ray.mu > 1.E-3 ? ray.mu : 1.E-3;
  ds_.bc.phi0 = ray.phi;
  if (ds_.flag.planck) {
    ds_.bc.btemp = ds_.temper[ds_.nlyr];
    ds_.bc.ttemp = 0.;
  }

  if (pmy_band_->TestFlag(RadiationFlags::BroadBand)) {
    // stellar source function overrides fbeam
    if (pmy_band_->TestFlag(RadiationFlags::StellarBeam)) {
      ds_.bc.fbeam = pmb->pimpl->planet->ParentInsolationFlux(wmin, wmax, 1.);
    } else if (pmy_band_->HasPar("fbeam_K")) {
      Real Ttop = pmy_band_->GetPar<Real>("fbeam_K");
      ds_.bc.fbeam = Constants::stefanBoltzmann * pow(Ttop, 4);
    } else if (pmy_band_->HasPar("fbeam")) {
      ds_.bc.fbeam = pmy_band_->GetPar<Real>("fbeam");
    } else {
      ds_.bc.fbeam = 0.;
    }
    ds_.bc.fbeam /= dist_au * dist_au;

    // planck source function
    ds_.wvnmlo = wmin;
    ds_.wvnmhi = wmax;
  }

  if (pmb != nullptr) {
    pmb->pcoord->Face1Area(k, j, pmb->is, pmb->ie + 1, farea_);
    pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol_);
  }

  //! \todo update this
  ds_.bc.fluor = fluor;
  ds_.bc.fisot = fisot;
  ds_.bc.albedo = albedo;
  ds_.bc.temis = 1.;

  auto &&uphi_umu = RadiationHelper::get_direction_grids(pmy_band_->rayOutput_);
  auto &uphi = uphi_umu.first;
  auto &umu = uphi_umu.second;

  if (umu.size() <= ds_.numu) {
    SetUserCosinePolarAngle(umu);

    for (int i = 0; i < umu.size(); ++i) {
      dir_axis_[i] = umu[i];
    }

    for (int i = umu.size(); i < ds_.numu; ++i) {
      dir_axis_[i] = 1.;
    }
  } else {
    throw RuntimeError("RTSolverDisort::Prepare",
                       "Number of polar angles in Disort is too small");
  }

  if (uphi.size() <= ds_.nphi) {
    SetUserAzimuthalAngle(uphi);

    for (int i = 0; i < uphi.size(); ++i) {
      dir_axis_[ds_.numu + i] = uphi[i];
    }

    for (int i = uphi.size(); i < ds_.nphi; ++i) {
      dir_axis_[ds_.numu + i] = 2. * M_PI;
    }
  } else {
    throw RuntimeError("RTSolverDisort::Prepare",
                       "Number of azimuthal angles in Disort is too small");
  }
}

void RadiationBand::RTSolverDisort::CalBandFlux(MeshBlock const *pmb, int k,
                                                int j, int il, int iu) {
  Real dist_au;

  if (pmy_band_->TestFlag(RadiationFlags::StellarBeam)) {
    Real time = pmb->pmy_mesh->time;
    if (pmy_band_->TestFlag(RadiationFlags::TimeDependent)) {
      dist_au = pmb->pimpl->planet->ParentDistanceInAu(time);
    } else {
      dist_au = pmb->pimpl->GetDistanceInAu();
    }
  }

  // bflxup and bflxdn has been reset in RadiationBand::CalBandFlux

  // loop over spectral grids in the band
  int b = 0;
  for (auto &spec : pmy_band_->pgrid_->spec) {
    if (!pmy_band_->TestFlag(RadiationFlags::BroadBand)) {
      // stellar source function
      if (pmy_band_->TestFlag(RadiationFlags::StellarBeam)) {
        ds_.bc.fbeam = pmb->pimpl->planet->ParentInsolationFlux(
            spec.wav1, spec.wav2, dist_au);
      } else {
        ds_.bc.fbeam = 0.;
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
    if (pmb != nullptr) {
      addDisortFlux(pmb->pcoord, b++, k, j, il, iu);
    } else {
      addDisortFlux(nullptr, b++, k, j, il, iu);
    }
  }
}

void RadiationBand::RTSolverDisort::addDisortFlux(Coordinates const *pcoord,
                                                  int b, int k, int j, int il,
                                                  int iu) {
  auto &bflxup = pmy_band_->bflxup;
  auto &bflxdn = pmy_band_->bflxdn;

  auto &flxup = pmy_band_->flxup_;
  auto &flxdn = pmy_band_->flxdn_;
  auto const &spec = pmy_band_->pgrid_->spec;

  /// accumulate flux from spectral bins
  for (int i = il; i <= iu; ++i) {
    int m = ds_.nlyr - (myrank_in_column_ * (iu - il) + i - il);
    //! \bug does not work for spherical geometry, need to scale area using
    //! farea(il)/farea(i)
    // flux up
    flxup(b, k, j, i) = ds_out_.rad[m].flup;

    //! \bug does not work for spherical geomtry, need to scale area using
    //! farea(il)/farea(i)
    // flux down
    flxdn(b, k, j, i) = ds_out_.rad[m].rfldir + ds_out_.rad[m].rfldn;

    bflxup(k, j, i) += spec[b].wght * flxup(b, k, j, i);
    bflxdn(k, j, i) += spec[b].wght * flxdn(b, k, j, i);
  }

  //! \note Spherical correction by XIZ
  //! xiz 2022 flux scaling so that the heating rate is the same as the
  //! plane-parallel scheme volheating scaling: first calculate flux divergence
  //! from DISORT using Plane-parallel in a cell then mulitpled by the cell
  //! volume divided by dx1f then solve for F using F1*S1-F2*S2 = volheating
  //! the top fluxes are the still the same as the plane-paralell values
  Real volh;
  Real bflxup_iu = bflxup(k, j, iu);
  Real bflxdn_iu = bflxdn(k, j, iu);

  if (pcoord != nullptr) {
    for (int i = iu - 1; i >= il; --i) {
      // upward
      volh = (bflxup_iu - bflxup(k, j, i)) / pcoord->dx1f(i) * vol_(i);
      bflxup_iu = bflxup(k, j, i);
      bflxup(k, j, i) =
          (bflxup(k, j, i + 1) * farea_(i + 1) - volh) / farea_(i);

      // downward
      volh = (bflxdn_iu - bflxdn(k, j, i)) / pcoord->dx1f(i) * vol_(i);
      bflxdn_iu = bflxdn(k, j, i);
      bflxdn(k, j, i) =
          (bflxdn(k, j, i + 1) * farea_(i + 1) - volh) / farea_(i);
    }
  }
}

void RadiationBand::RTSolverDisort::CalBandRadiance(MeshBlock const *pmb, int k,
                                                    int j) {
  if (ds_.flag.onlyfl) {
    throw RuntimeError("RTSolverDisort::CalBandRadiance",
                       "Radiance calculation disabled");
  }

  if (ds_.ntau != 1) {
    throw RuntimeError("RTSolverDisort::CalBandRadiance",
                       "Only toa radiance (ds.ntau = 1) is supported");
  }

  int nrays = ds_.nphi * ds_.numu;

  if (nrays < pmy_band_->GetNumOutgoingRays()) {
    throw RuntimeError("RTSolverDisort::CalBandRadiance",
                       "Number of outgoing rays more than DISORT can host");
  }

  // toa has been reset in RadiationBand::CalBandRadiance

  Real dist_au;
  if (pmy_band_->TestFlag(RadiationFlags::StellarBeam)) {
    if (pmy_band_->TestFlag(RadiationFlags::TimeDependent)) {
      Real time = pmb->pmy_mesh->time;
      dist_au = pmb->pimpl->planet->ParentDistanceInAu(time);
    } else {
      dist_au = pmb->pimpl->GetDistanceInAu();
    }
  }

  RadiationHelper::get_phase_momentum(ds_.pmom, ISOTROPIC, 0., ds_.nmom);

  // loop over spectral grids in the band
  int b = 0;
  for (auto &spec : pmy_band_->pgrid_->spec) {
    // override source function for non-broadband radiation
    if (!(pmy_band_->TestFlag(RadiationFlags::BroadBand))) {
      // stellar source function
      if (pmy_band_->TestFlag(RadiationFlags::StellarBeam)) {
        ds_.bc.fbeam = pmb->pimpl->planet->ParentInsolationFlux(
            spec.wav1, spec.wav2, dist_au);
      } else {
        ds_.bc.fbeam = 0.0;
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
    addDisortRadiance(b++, k, j);
  }
}

void RadiationBand::RTSolverDisort::addDisortRadiance(int b, int k, int j) {
  auto &toa = pmy_band_->toa_;
  auto &btoa = pmy_band_->btoa;
  auto &spec = pmy_band_->pgrid_->spec;
  auto &rayout = pmy_band_->rayOutput_;

  for (int n = 0; n < pmy_band_->GetNumOutgoingRays(); ++n) {
    Real val;
    Real coor[2] = {rayout[n].mu, rayout[n].phi};
    interpn(&val, coor, ds_out_.uu, dir_axis_.data(), dir_dim_, 2, 1);
    toa(b, n, k, j) = val;
    btoa(n, k, j) += spec[b].wght * toa(b, n, k, j);
  }
}

#endif  // RT_DISORT
