// C/C++ headers
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

// external
#include <yaml-cpp/yaml.h>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// climath
#include <climath/core.h>

// utils
#include <utils/vectorize.hpp>

// astro
#include <astro/celestrial_body.hpp>

// canoe
#include <common.hpp>
#include <impl.hpp>

// harp
#include "radiation.hpp"
#include "radiation_band.hpp"
#include "rt_solvers.hpp"

const std::string Radiation::input_key = "radiation_config";

Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin) {
  Application::Logger app("harp");
  app->Log("Initialize Radiation");

  // radiation flags
  SetFlag(RadiationHelper::parse_radiation_flags(
      pin->GetOrAddString("radiation", "flags", "")));

  // radiation bands
  bands_ = RadiationBandsFactory::CreateFrom(pin, input_key);

  // incoming radiation direction (mu, phi) in degree
  std::string str = pin->GetOrAddString("radiation", "indir", "(0.,0.)");
  rayInput_ = RadiationHelper::parse_radiation_directions(str);

  int ncells1 = pmb->ncells1;
  int ncells2 = pmb->ncells2;
  int ncells3 = pmb->ncells3;

  // output radiance
  int nout = GetNumOutgoingRays();
  if (nout > 0) {
    radiance.NewAthenaArray(nout, ncells3, ncells2);
    // set band toa
    int n = 0;
    for (auto &p : bands_) {
      p->btoa.InitWithShallowSlice(radiance, 3, n, p->GetNumOutgoingRays());
      n += p->GetNumOutgoingRays();
    }
  }

  // output flux
  flxup.NewAthenaArray(bands_.size(), ncells3, ncells2, ncells1 + 1);
  flxdn.NewAthenaArray(bands_.size(), ncells3, ncells2, ncells1 + 1);

  for (int n = 0; n < bands_.size(); ++n) {
    bands_[n]->bflxup.InitWithShallowSlice(flxup, 4, n, 1);
    bands_[n]->bflxup.InitWithShallowSlice(flxdn, 4, n, 1);
  }

  // time control
  SetCooldownTime(pin->GetOrAddReal("radiation", "dt", 0.));
}

Radiation::~Radiation() {
  Application::Logger app("harp");
  app->Log("Destroy Radiation");
}

RadiationBandPtr Radiation::GetBandByName(std::string const &name) const {
  for (auto &band : bands_) {
    if (band->GetName() == name) {
      return band;
    }
  }
  throw NotFoundError("GetBand", "Band " + name);
}

size_t Radiation::GetNumOutgoingRays() const {
  size_t num = 0;
  for (auto &p : bands_) {
    num += p->GetNumOutgoingRays();
  }
  return num;
}

void Radiation::CalRadiativeFlux(MeshBlock const *pmb, Real time, int k, int j,
                                 int il, int iu) {
  auto pcoord = pmb->pcoord;
  auto planet = pmb->pimpl->planet;
  Real dist = pmb->pimpl->GetDistanceInAu();

  Direction ray;
  if (TestFlag(RadiationFlags::Dynamic)) {
    ray = planet->ParentZenithAngle(time, pcoord->x2v(j), pcoord->x3v(k));
    dist = planet->ParentDistanceInAu(time);
  } else {
    ray = rayInput_[0];
  }

  AirColumn &&ac =
      AirParcelHelper::gather_from_primitive(pmb, k, j, 0, pmb->ncells1 - 1);
  for (auto &air : ac) air.ToMoleFraction();

  for (auto &p : bands_) {
    // iu ~= ie + 1
    p->SetSpectralProperties(ac, pcoord, k, j);
    p->psolver->CalBandFlux(ray, dist, k, j, il, iu);
  }
}

void Radiation::CalRadiance(MeshBlock const *pmb, Real time, int k, int j,
                            int il, int iu) {
  Application::Logger app("harp");
  app->Log("CalRadiance");

  auto pcoord = pmb->pcoord;
  auto planet = pmb->pimpl->planet;
  Real dist = pmb->pimpl->GetDistanceInAu();

  Direction ray;
  if (TestFlag(RadiationFlags::Dynamic)) {
    ray = planet->ParentZenithAngle(time, pcoord->x2v(j), pcoord->x3v(k));
    dist = planet->ParentDistanceInAu(time);
  } else {
    ray = rayInput_[0];
  }

  AirColumn &&ac =
      AirParcelHelper::gather_from_primitive(pmb, k, j, 0, pmb->ncells1 - 1);
  for (auto &air : ac) air.ToMoleFraction();

  for (auto &p : bands_) {
    // iu ~= ie + 1
    p->SetSpectralProperties(ac, pcoord, k, j);
    p->psolver->CalBandRadiance(ray, dist, k, j, il, iu);
  }
}

void Radiation::AddRadiativeFlux(Hydro *phydro, int k, int j, int il,
                                 int iu) const {
  // x1-flux divergence
  for (size_t b = 0; b < bands_.size(); ++b) {
#pragma omp simd
    for (int i = il; i <= iu; ++i)
      phydro->flux[X1DIR](IEN, k, j, i) +=
          flxup(b, k, j, i) - flxdn(b, k, j, i);
  }
}

size_t Radiation::RestartDataSizeInBytes() const {
  return flxup.GetSizeInBytes() + flxdn.GetSizeInBytes();
}

size_t Radiation::DumpRestartData(char *pdst) const {
  int offset = 0;

  std::memcpy(pdst + offset, flxup.data(), flxup.GetSizeInBytes());
  offset += flxup.GetSizeInBytes();
  std::memcpy(pdst + offset, flxdn.data(), flxdn.GetSizeInBytes());
  offset += flxdn.GetSizeInBytes();

  return RestartDataSizeInBytes();
}

size_t Radiation::LoadRestartData(char *psrc) {
  Application::Logger app("harp");
  app->Log("Radiation restarted");
  int offset = 0;

  std::memcpy(flxup.data(), psrc + offset, flxup.GetSizeInBytes());
  offset += flxup.GetSizeInBytes();
  std::memcpy(flxdn.data(), psrc + offset, flxdn.GetSizeInBytes());
  offset += flxdn.GetSizeInBytes();

  return RestartDataSizeInBytes();
}

namespace RadiationHelper {
std::vector<Direction> parse_radiation_directions(std::string str) {
  std::vector<std::string> dstr = Vectorize<std::string>(str.c_str());
  int nray = dstr.size();

  std::vector<Direction> ray(nray);

  auto jt = dstr.begin();
  for (auto it = ray.begin(); it != ray.end(); ++it, ++jt) {
    it->phi = 0.;
    sscanf(jt->c_str(), "(%lf,%lf)", &it->mu, &it->phi);
    it->mu = cos(deg2rad(it->mu));
    it->phi = deg2rad(it->phi);
  }

  return ray;
}

uint64_t parse_radiation_flags(std::string str) {
  uint64_t flags = 0LL;
  std::stringstream msg;

  std::vector<std::string> dstr = Vectorize<std::string>(str.c_str(), " ,");
  for (int i = 0; i < dstr.size(); ++i) {
    if (dstr[i] == "static") {
      flags &= !RadiationFlags::Dynamic;
    } else if (dstr[i] == "dynamic") {
      flags |= RadiationFlags::Dynamic;
    } else if (dstr[i] == "bin") {
      flags &= !RadiationFlags::LineByLine;
    } else if (dstr[i] == "lbl") {
      flags |= RadiationFlags::LineByLine;
    } else if (dstr[i] == "ck") {
      flags |= RadiationFlags::CorrelatedK;
    } else if (dstr[i] == "planck") {
      flags |= RadiationFlags::Planck;
    } else if (dstr[i] == "star") {
      flags |= RadiationFlags::Star;
    } else if (dstr[i] == "spher") {
      flags |= RadiationFlags::Sphere;
    } else if (dstr[i] == "only") {
      flags |= RadiationFlags::FluxOnly;
    } else if (dstr[i] == "normalize") {
      flags |= RadiationFlags::Normalize;
    } else if (dstr[i] == "write_bin_radiance") {
      flags |= RadiationFlags::WriteBinRadiance;
    } else {
      msg << "flag:" << dstr[i] << "unrecognized" << std::endl;
      throw RuntimeError("parse_radiation_flags", msg.str());
    }

    // check flags consistency
    if ((flags & RadiationFlags::LineByLine) &&
        (flags & RadiationFlags::CorrelatedK)) {
      msg << "ck cannot be used with lbl." << std::endl;
      throw RuntimeError("parse_radiation_flags", msg.str());
    }
  }

  return flags;
}

void get_phase_momentum(Real *pmom, int iphas, Real gg, int npmom) {
  pmom[0] = 1.;
  for (int k = 1; k < npmom; k++) pmom[k] = 0.;

  switch (iphas) {
    case 0:  // HENYEY_GREENSTEIN
      if (gg <= -1. || gg >= 1.)
        std::cout << "getmom--bad input variable gg" << std::endl;
      for (int k = 1; k <= npmom; k++) pmom[k] = pow(gg, (Real)k);
      break;

    case 1:  // RAYLEIGH
      pmom[2] = 0.1;
      break;
  }
}

}  // namespace RadiationHelper
