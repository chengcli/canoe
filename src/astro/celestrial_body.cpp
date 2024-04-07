// C/C++
#include <sstream>
#include <stdexcept>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// climath
#include <climath/core.h>
#include <climath/interpolation.h>

// harp
#include <harp/spectral_grid.hpp>

// utils
#include <utils/fileio.hpp>

// astro
#include "celestrial_body.hpp"

void CelestrialBody::loadOrbitalData(ParameterInput *pin) {
  char entry[80];
  auto name = GetName();
  snprintf(entry, sizeof(entry), "%s.re", name.c_str());
  re = km2m(pin->GetOrAddReal("astronomy", entry, 0.));

  snprintf(entry, sizeof(entry), "%s.rp", name.c_str());
  rp = km2m(pin->GetOrAddReal("astronomy", entry, re));

  snprintf(entry, sizeof(entry), "%s.obliq", name.c_str());
  obliq = deg2rad(pin->GetOrAddReal("astronomy", entry, 0.));

  snprintf(entry, sizeof(entry), "%s.spinp", name.c_str());
  spinp = day2sec(pin->GetOrAddReal("astronomy", entry, 0.));

  snprintf(entry, sizeof(entry), "%s.orbit_a", name.c_str());
  orbit_a = au2m(pin->GetOrAddReal("astronomy", entry, 0.));

  snprintf(entry, sizeof(entry), "%s.orbit_e", name.c_str());
  orbit_e = pin->GetOrAddReal("astronomy", entry, 0.);

  snprintf(entry, sizeof(entry), "%s.orbit_i", name.c_str());
  orbit_i = deg2rad(pin->GetOrAddReal("astronomy", entry, 0.));

  snprintf(entry, sizeof(entry), "%s.orbit_p", name.c_str());
  orbit_p = day2sec(pin->GetOrAddReal("astronomy", entry, 0.));

  snprintf(entry, sizeof(entry), "%s.equinox", name.c_str());
  equinox = pin->GetOrAddReal("astronomy", entry, 0.);

  snprintf(entry, sizeof(entry), "%s.grav_eq", name.c_str());
  grav_eq = pin->GetOrAddReal("astronomy", entry, 0.);
}

CelestrialBody::CelestrialBody(ParameterInput *pin, std::string name)
    : NamedGroup(name) {
  Application::Logger app("astro");
  app->Log("Initialize CelestrialBody " + GetName());

  loadOrbitalData(pin);

  if (pin->DoesParameterExist("astronomy", name + ".parent")) {
    std::string parent_name = pin->GetString("astronomy", name + ".parent");
    parent_ = std::make_shared<CelestrialBody>(pin, parent_name);
  }

  if (pin->DoesParameterExist("astronomy", name + ".spectra")) {
    std::string sfile = pin->GetString("astronomy", name + ".spectra");
    if (!FileExists(sfile)) {
      app->Error("Cannot open spectral file " + sfile);
    } else {
      loadSpectralData(sfile);
    }
  }
}

CelestrialBody::~CelestrialBody() {
  Application::Logger app("astro");
  app->Log("Destroy CelestrialBody" + GetName());
}

void CelestrialBody::loadSpectralData(std::string sfile) {
  AthenaArray<Real> spectrum;
  ReadDataTable(&spectrum, sfile);

  int nspec = spectrum.GetDim2();
  spec_.resize(nspec);

  for (int i = 0; i < nspec; ++i) {
    spec_[i].x = spectrum(i, 0);
    spec_[i].y = spectrum(i, 1);
  }

  spline(spec_.size(), spec_.data(), 0., 0.);
}

Direction CelestrialBody::ParentZenithAngle(Real time, Real lat,
                                            Real lon) const {
  Direction dir;

  if (spinp == 0. && orbit_p == 0.) dir.mu = cos(-lon + M_PI) * cos(lat);
  if (spinp == 0. && orbit_p != 0.)
    dir.mu = cos((time * 2. * M_PI * (-1. / orbit_p)) - lon + M_PI) * cos(lat);
  if (spinp != 0. && orbit_p == 0.)
    dir.mu = cos((time * 2. * M_PI * (1. / spinp)) - lon + M_PI) * cos(lat);
  if (spinp != 0. && orbit_p != 0.)
    dir.mu =
        cos((time * 2. * M_PI * (1. / spinp - 1. / orbit_p)) - lon + M_PI) *
        cos(lat);
  dir.phi = 0.;

  return dir;
}

Real CelestrialBody::ParentInsolationFlux(Real wav, Real dist_au) const {
  if (!parent_) {
    throw RuntimeError("ParentInsolationFlux", "no parent body");
  }

  Real dx;
  // assuming ascending in wav
  if ((il_ >= 0) && (wav < parent_->spec_[il_].x)) il_ = -1;
  il_ = find_place_in_table(parent_->spec_.size(), parent_->spec_.data(), wav,
                            &dx, il_);
  Real flux = splint(wav, parent_->spec_.data() + il_, dx);
  return flux / (dist_au * dist_au);
}

Real CelestrialBody::ParentInsolationFlux(Real wlo, Real whi,
                                          Real dist_au) const {
  if (!parent_) {
    throw RuntimeError("ParentInsolationFlux", "no parent body");
  }

  //! \todo check whether this is correct
  if (wlo == whi) return ParentInsolationFlux(wlo, dist_au);
  // assuming ascending in wav
  Real dx;
  int il = find_place_in_table(parent_->spec_.size(), parent_->spec_.data(),
                               wlo, &dx, -1);
  int iu = find_place_in_table(parent_->spec_.size(), parent_->spec_.data(),
                               whi, &dx, il);
  Real flux = 0.5 * parent_->spec_[il].y;
  for (int i = il; i < iu; ++i) {
    flux += 0.5 * (parent_->spec_[i].y + parent_->spec_[i + 1].y) *
            (parent_->spec_[i + 1].x - parent_->spec_[i].x);
  }
  return flux / (dist_au * dist_au);
}

Real CelestrialBody::ParentDistanceInAu(Real time) const {
  Real orbit_b = sqrt(1. - orbit_e * orbit_e) * orbit_a;
  Real orbit_c = orbit_b * orbit_b / orbit_a;
  return m2au(orbit_c / (1. + orbit_e * cos(2. * M_PI / orbit_p * time -
                                            equinox - M_PI / 2.)));
}

CelestrialBodyPtr PlanetFactory::CreateFrom(MeshBlock *pmb,
                                            ParameterInput *pin) {
  if (pin->DoesParameterExist("astronomy", "planet")) {
    std::string name = pin->GetString("astronomy", "planet");
    return std::make_shared<CelestrialBody>(pin, name);
  } else {
    return nullptr;
  }
}
