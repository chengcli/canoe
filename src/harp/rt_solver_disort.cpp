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

// exo3
#include <exo3/cubed_sphere.hpp>

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
    // SetFlags is a member function of DisortWrapper
    SetFlags(to_map_bool(rad["Disort-flags"]));
  }

  if (pmy_band->HasPar("btemp")) {
    ds_.bc.btemp = pmy_band->GetPar<Real>("btemp");
  }

  if (pmy_band->HasPar("albedo")) {
    ds_.bc.albedo = pmy_band->GetPar<Real>("albedo");
  } else {
    ds_.bc.albedo = 1.;
  }

  if (pmy_band->HasPar("ttemp")) {
    ds_.bc.ttemp = pmy_band->GetPar<Real>("ttemp");
  }

  if (pmy_band->HasPar("temis")) {
    ds_.bc.temis = pmy_band->GetPar<Real>("temis");
  } else {
    ds_.bc.temis = 0.;
  }

  if (pmy_band->HasPar("fluor")) {
    ds_.bc.fluor = pmy_band->GetPar<Real>("fluor");
  } else {
    ds_.bc.fluor = 0.;
  }

  if (pmy_band->HasPar("fisot")) {
    ds_.bc.fisot = pmy_band->GetPar<Real>("fisot");
  } else {
    ds_.bc.fisot = 0.;
  }

  Application::Logger app("harp");

  // override disort planck flag
  if (pmy_band->TestFlag(RadiationFlags::ThermalEmission)) {
    ds_.flag.planck = true;
    app->Log("Planck function is enabled for band " + pmy_band->GetName());
  } else {
    ds_.flag.planck = false;
    app->Log("Planck function is disabled for band " + pmy_band->GetName());
  }

  SetHeader("Disort solving band " + pmy_band_->GetName());
}

//! \todo update based on band outdir
void RadiationBand::RTSolverDisort::Resize(int nlyr, int nstr) {
  RadiationBand::RTSolver::Resize(nlyr, nstr);
  Unseal();

  auto &rayout = pmy_band_->rayOutput_;
  auto &&uphi_umu = RadiationHelper::get_direction_grids(rayout);

  SetAtmosphereDimension(nlyr, nstr, nstr);

  dir_dim_[0] = uphi_umu.second.size();  // umu
  dir_dim_[1] = uphi_umu.first.size();   // uphi
  dir_axis_.resize(dir_dim_[0] + dir_dim_[1]);

  SetIntensityDimension(std::max(dir_dim_[1], 1lu), 1,
                        std::max(dir_dim_[0], 1lu));
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

  Real dist_au = 1.;
  Direction ray = pmb->pimpl->prad->GetRayInput(0);
  auto planet = pmb->pimpl->planet;

  if (planet && pmy_band_->TestFlag(RadiationFlags::TimeDependent)) {
    Real time = pmb->pmy_mesh->time;
    Real lat, lon, colat;
#ifdef CUBED_SPHERE
    pmb->pimpl->pexo3->GetLatLon(&lat, &lon, k, j, pmb->ie);
    colat = M_PI / 2. - lat;
#else   // FIXME: add another condition
    colat = pmb->pcoord->x2v(j);
    lon = pmb->pcoord->x3v(k);
#endif  // CUBED_SPHERE
    ray = planet->ParentZenithAngle(time, colat, lon);
    dist_au = planet->ParentDistanceInAu(time);
  } else {  // constant zenith angle
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

  pmy_band_->Regroup(pmb, X1DIR);
  myrank_in_column_ = pmy_band_->GetRankInGroup();

  // transfer temperature
  if (ds_.flag.planck) {
    pmy_band_->PackTemperature();
    pmy_band_->Transfer(pmb, 0);
    pmy_band_->UnpackTemperature(&ds_);
  }

  ds_.bc.umu0 = ray.mu > 1.E-3 ? ray.mu : 1.E-3;
  ds_.bc.phi0 = ray.phi;

  ds_.wvnmlo = wmin;
  ds_.wvnmhi = wmax;

  if (pmy_band_->TestFlag(RadiationFlags::BroadBand)) {
    // stellar source function overrides fbeam
    if (pmy_band_->HasPar("S0")) {
      ds_.bc.fbeam = pmy_band_->GetPar<Real>("S0");
    } else if (pmy_band_->HasPar("temp0")) {
      Real temp0 = pmy_band_->GetPar<Real>("temp0");
      ds_.bc.fbeam = Constants::stefanBoltzmann * pow(temp0, 4);
    } else if (planet && planet->HasParentFlux()) {
      ds_.bc.fbeam = planet->ParentInsolationFlux(wmin, wmax, 1.);
    } else {
      ds_.bc.fbeam = 0.;
    }
    ds_.bc.fbeam /= dist_au * dist_au;
  }

  if (pmb->pcoord != nullptr) {
    pmb->pcoord->Face1Area(k, j, pmb->is, pmb->ie + 1, farea_);
    pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol_);
  }

  auto &&uphi_umu = RadiationHelper::get_direction_grids(pmy_band_->rayOutput_);
  auto &uphi = uphi_umu.first;
  auto &umu = uphi_umu.second;

  if (umu.size() > 0) {
    SetUserCosinePolarAngle(umu);

    for (int i = 0; i < umu.size(); ++i) {
      dir_axis_[i] = umu[i];
    }
  }

  if (uphi.size() > 0) {
    SetUserAzimuthalAngle(uphi);

    for (int i = 0; i < uphi.size(); ++i) {
      dir_axis_[ds_.numu + i] = uphi[i];
    }
  }
}

void RadiationBand::RTSolverDisort::CalBandFlux(MeshBlock const *pmb, int k,
                                                int j, int il, int iu) {
  Real dist_au = 1.;
  auto planet = pmb->pimpl->planet;

  if (planet && pmy_band_->TestFlag(RadiationFlags::TimeDependent)) {
    dist_au = planet->ParentDistanceInAu(pmb->pmy_mesh->time);
  } else if (pmy_band_->HasPar("dist_au")) {
    dist_au = pmy_band_->GetPar<Real>("dist_au");
  }

  // bflxup and bflxdn has been reset in RadiationBand::CalBandFlux

  // loop over spectral grids in the band
  int b = 0;
  bool override_with_stellar_spectra = false;
  if (!pmy_band_->TestFlag(RadiationFlags::BroadBand) &&
      !pmy_band_->HasPar("S0") && !pmy_band_->HasPar("temp0") && planet &&
      planet->HasParentFlux()) {
    override_with_stellar_spectra = true;
  }

  for (auto &spec : pmy_band_->pgrid_->spec) {
    if (override_with_stellar_spectra) {
      // stellar source function
      ds_.bc.fbeam =
          planet->ParentInsolationFlux(spec.wav1, spec.wav2, dist_au);
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
    if (pmb->pcoord != nullptr) {
      addDisortFlux(pmb->pcoord, b++, k, j, il, iu);
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

  Real dist_au = 1.;
  auto planet = pmb->pimpl->planet;

  if (planet && pmy_band_->TestFlag(RadiationFlags::TimeDependent)) {
    dist_au = pmb->pimpl->planet->ParentDistanceInAu(pmb->pmy_mesh->time);
  } else if (pmy_band_->HasPar("dist_au")) {
    dist_au = pmy_band_->GetPar<Real>("dist_au");
  }

  RadiationHelper::get_phase_momentum(ds_.pmom, ISOTROPIC, 0., ds_.nmom);

  // loop over spectral grids in the band
  int b = 0;
  for (auto &spec : pmy_band_->pgrid_->spec) {
    // override source function for non-broadband radiation
    if (!(pmy_band_->TestFlag(RadiationFlags::BroadBand))) {
      // stellar source function
      if (planet) {
        ds_.bc.fbeam =
            planet->ParentInsolationFlux(spec.wav1, spec.wav2, dist_au);
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
