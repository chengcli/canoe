// C/C++ headers
#include <climath/core.h>
#include <climath/interpolation.h>

#include <athena/athena.hpp>
#include <athena/parameter_input.hpp>
#include <debugger/debugger.hpp>
#include <sstream>
#include <stdexcept>
#include <utils/fileio.hpp>

// cliastro headers
#include "celestrial_body.hpp"
// #include "../radiation/radiation.hpp"

void CelestrialBody::ReadCelestrialData_(ParameterInput *pin,
                                         std::string myname) {
  char entry[80];
  snprintf(entry, 80, "%s.re", name.c_str());
  re = km2m(pin->GetOrAddReal("astronomy", entry, 0.));

  snprintf(entry, 80, "%s.rp", name.c_str());
  rp = km2m(pin->GetOrAddReal("astronomy", entry, re));

  snprintf(entry, 80, "%s.obliq", name.c_str());
  obliq = deg2rad(pin->GetOrAddReal("astronomy", entry, 0.));

  snprintf(entry, 80, "%s.spinp", name.c_str());
  spinp = day2sec(pin->GetOrAddReal("astronomy", entry, 0.));

  snprintf(entry, 80, "%s.orbit_a", name.c_str());
  orbit_a = au2m(pin->GetOrAddReal("astronomy", entry, 0.));

  snprintf(entry, 80, "%s.orbit_e", name.c_str());
  orbit_e = pin->GetOrAddReal("astronomy", entry, 0.);

  snprintf(entry, 80, "%s.orbit_i", name.c_str());
  orbit_i = deg2rad(pin->GetOrAddReal("astronomy", entry, 0.));

  snprintf(entry, 80, "%s.orbit_p", name.c_str());
  orbit_p = day2sec(pin->GetOrAddReal("astronomy", entry, 0.));

  snprintf(entry, 80, "%s.equinox", name.c_str());
  equinox = pin->GetOrAddReal("astronomy", entry, 0.);

  snprintf(entry, 80, "%s.grav_eq", name.c_str());
  grav_eq = pin->GetOrAddReal("astronomy", entry, 0.);
}

CelestrialBody::CelestrialBody(ParameterInput *pin)
    : parent(nullptr), spec_(nullptr), il_(-1) {
  std::stringstream msg;
  name = pin->GetOrAddString("astronomy", "planet", "unknown");
  ReadCelestrialData_(pin, name);

  if (pin->DoesParameterExist("astronomy", name + ".parent")) {
    std::string parent_name = pin->GetString("astronomy", name + ".parent");
    parent = new CelestrialBody(pin, parent_name);
  } else
    parent = nullptr;

  if (pin->DoesParameterExist("astronomy", name + ".spec_file")) {
    std::string sfile = pin->GetString("astronomy", name + ".spec_file");
    if (!FileExists(sfile)) {
      Debugger::Fatal("CelestrialBody", "Cannot open spectral file", sfile);
    } else {
      ReadSpectraFile(sfile);
    }
  } else {
    spec_ = new float_triplet[1];
    nspec_ = 1;
  }
}

CelestrialBody::CelestrialBody(ParameterInput *pin, std::string myname)
    : parent(nullptr), name(myname), spec_(nullptr), il_(-1) {
  std::stringstream msg;
  ReadCelestrialData_(pin, name);

  if (pin->DoesParameterExist("astronomy", name + ".parent")) {
    std::string parent_name = pin->GetString("astronomy", name + ".parent");
    parent = new CelestrialBody(pin, parent_name);
  } else
    parent = nullptr;

  if (pin->DoesParameterExist("astronomy", name + ".spec_file")) {
    std::string sfile = pin->GetString("astronomy", name + ".spec_file");
    if (!FileExists(sfile)) {
      Debugger::Fatal("CelestrialBody", "Cannot open spectral file", sfile);
    } else {
      ReadSpectraFile(sfile);
    }
  } else {
    spec_ = new float_triplet[1];
    nspec_ = 1;
  }
}

CelestrialBody::~CelestrialBody() {
  if (parent != nullptr) delete parent;
  delete[] spec_;
}

void CelestrialBody::ReadSpectraFile(std::string sfile) {
  AthenaArray<Real> spectrum;
  ReadDataTable(spectrum, sfile);

  nspec_ = spectrum.GetDim2();

  if (spec_ != nullptr) delete[] spec_;
  spec_ = new float_triplet[nspec_];
  for (int i = 0; i < nspec_; ++i) {
    spec_[i].x = spectrum(i, 0);
    spec_[i].y = spectrum(i, 1);
  }

  spline(nspec_, spec_, 0., 0.);
}

void CelestrialBody::ParentZenithAngle(Real *mu, Real *phi, Real time,
                                       Real colat, Real lon) {
  Real lat = M_PI / 2. - colat;
  if (spinp == 0. && orbit_p == 0.) *mu = cos(-lon + M_PI) * cos(lat);
  if (spinp == 0. && orbit_p != 0.)
    *mu = cos((time * 2. * M_PI * (-1. / orbit_p)) - lon + M_PI) * cos(lat);
  if (spinp != 0. && orbit_p == 0.)
    *mu = cos((time * 2. * M_PI * (1. / spinp)) - lon + M_PI) * cos(lat);
  if (spinp != 0. && orbit_p != 0.)
    *mu = cos((time * 2. * M_PI * (1. / spinp - 1. / orbit_p)) - lon + M_PI) *
          cos(lat);
  *phi = 0.;
}

Real CelestrialBody::ParentInsolationFlux(Real wav, Real dist_au) {
  Real dx;
  // assuming ascending in wav
  if ((il_ >= 0) && (wav < parent->spec_[il_].x)) il_ = -1;
  il_ = find_place_in_table(parent->nspec_, parent->spec_, wav, &dx, il_);
  Real flux = splint(wav, parent->spec_ + il_, dx);
  return flux / (dist_au * dist_au);
}

Real CelestrialBody::ParentInsolationFlux(Real wlo, Real whi, Real dist_au) {
  //! \todo check whether this is correct
  if (wlo == whi) return ParentInsolationFlux(wlo, dist_au);
  // assuming ascending in wav
  Real dx;
  int il = find_place_in_table(parent->nspec_, parent->spec_, wlo, &dx, -1);
  int iu = find_place_in_table(parent->nspec_, parent->spec_, whi, &dx, il);
  Real flux = 0.5 * parent->spec_[il].y;
  for (int i = il; i < iu; ++i) {
    flux += 0.5 * (parent->spec_[i].y + parent->spec_[i + 1].y) *
            (parent->spec_[i + 1].x - parent->spec_[i].x);
  }
  return flux / (dist_au * dist_au);
}

Real CelestrialBody::ParentDistanceInAu(Real time) {
  Real orbit_b = sqrt(1. - orbit_e * orbit_e) * orbit_a;
  Real orbit_c = orbit_b * orbit_b / orbit_a;
  return m2au(orbit_c / (1. + orbit_e * cos(2. * M_PI / orbit_p * time -
                                            equinox - M_PI / 2.)));
}
