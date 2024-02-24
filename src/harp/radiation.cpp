// C/C++ headers
#include <algorithm>
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

  // dimensions
  int ncells1 = pmb->ncells1;
  int ncells2 = pmb->ncells2;
  int ncells3 = pmb->ncells3;

  // radiation flags
  SetFlag(RadiationHelper::parse_radiation_flags(
      pin->GetOrAddString("radiation", "flags", "")));

  // radiation bands
  bands_ = RadiationBandsFactory::CreateFrom(pin, input_key);

  // incoming radiation direction (mu, phi) in degrees
  std::string str = pin->GetOrAddString("radiation", "indir", "(0.,0.)");
  rayInput_ = RadiationHelper::parse_radiation_directions(str);

  // radiation configuration
  int nstr = pin->GetOrAddInteger("radiation", "nstr", 8);

  for (auto &p : bands_) {
    // outgoing radiation direction (mu,phi) in degrees
    if (pin->DoesParameterExist("radiation", p->GetName() + ".outdir")) {
      auto str = pin->GetString("radiation", p->GetName() + ".outdir");
      p->SetOutgoingRays(RadiationHelper::parse_radiation_directions(str));
    } else if (pin->DoesParameterExist("radiation", "outdir")) {
      auto str = pin->GetString("radiation", "outdir");
      p->SetOutgoingRays(RadiationHelper::parse_radiation_directions(str));
    }

    // allocate memory
    p->Resize(ncells1, ncells2, ncells3, nstr);
  }

  // output radiance
  int nout = GetNumOutgoingRays();
  if (nout > 0) {
    radiance.NewAthenaArray(nout, ncells3, ncells2);
    // set band toa
    int n = 0;
    for (auto &p : bands_) {
      //! Delete the old array and initialize with a shallow slice
      p->btoa.DeleteAthenaArray();
      p->btoa.InitWithShallowSlice(radiance, 3, n, p->GetNumOutgoingRays());
      n += p->GetNumOutgoingRays();
    }
  }

  // output flux
  flxup.NewAthenaArray(bands_.size(), ncells3, ncells2, ncells1 + 1);
  flxdn.NewAthenaArray(bands_.size(), ncells3, ncells2, ncells1 + 1);

  for (int n = 0; n < bands_.size(); ++n) {
    //! Delete the old array and initialize with a shallow slice
    bands_[n]->bflxup.DeleteAthenaArray();
    bands_[n]->bflxup.InitWithShallowSlice(flxup, 4, n, 1);

    //! Delete the old array and initialize with a shallow slice
    bands_[n]->bflxdn.DeleteAthenaArray();
    bands_[n]->bflxdn.InitWithShallowSlice(flxdn, 4, n, 1);
  }

  // radiative time scale
  rtime.NewAthenaArray(ncells3, ncells2, ncells1);

  // time control
  SetCooldownTime(pin->GetOrAddReal("radiation", "dt", 0.));

  // relaxation time
  relax_time_ = pin->GetOrAddReal("radiation", "relax_time", 1.);
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

void Radiation::CalFlux(MeshBlock const *pmb, int k, int j, int il, int iu) {
  auto pcoord = pmb->pcoord;

  AirColumn &&ac = AirParcelHelper::gather_from_primitive(pmb, k, j);

  Real grav = -pmb->phydro->hsrc.GetG1();
  Real H0 = pmb->pcoord->GetPressureScaleHeight();

  for (auto &p : bands_) {
    // iu ~= ie + 1
    p->SetSpectralProperties(ac, pcoord->x1f.data(), grav * H0, k, j);
    p->CalBandFlux(pmb, k, j, il, iu);
  }
}

void Radiation::CalRadiance(MeshBlock const *pmb, int k, int j) {
  auto pcoord = pmb->pcoord;

  AirColumn &&ac = AirParcelHelper::gather_from_primitive(pmb, k, j);

  Real grav = -pmb->phydro->hsrc.GetG1();
  Real H0 = pmb->pcoord->GetPressureScaleHeight();

  for (auto &p : bands_) {
    // iu ~= ie + 1
    p->SetSpectralProperties(ac, pcoord->x1f.data(), grav * H0, k, j);
    p->CalBandRadiance(pmb, k, j);
  }
}

void Radiation::CalTimeStep(MeshBlock const* pmb, int k, int j, int il,
                            int iu) {
  Real total_flux1 = 0., total_flux2 = 0.;
  auto pcoord = pmb->pcoord;
  auto phydro = pmb->phydro;

  time_step_ = 1.e99;

  for (size_t b = 0; b < bands_.size(); ++b) {
    total_flux1 += flxup(b, k, j, il) - flxdn(b, k, j, il);
  }

  for (int i = il; i <= iu; ++i) {
    for (size_t b = 0; b < bands_.size(); ++b) {
      total_flux2 += flxup(b, k, j, i + 1) - flxdn(b, k, j, i + 1);
    }
    rtime(k, j, i) = pcoord->dx1f(i) * phydro->u(IEN, k, j, i) 
      / (total_flux2 - total_flux1);
    total_flux1 = total_flux2;
    if (rtime(k, j, i) > 0.) {
      time_step_ = std::min(time_step_, rtime(k, j, i));
    }
  }
}

void Radiation::AddRadiativeFlux(Hydro *phydro, int k, int j, int il,
                                 int iu) const {
  // x1-flux divergence
  for (size_t b = 0; b < bands_.size(); ++b) {
    for (int i = il; i <= iu; ++i)
      phydro->flux[X1DIR](IEN, k, j, i) +=
          flxup(b, k, j, i) - flxdn(b, k, j, i);
  }
}

size_t Radiation::RestartDataSizeInBytes(Mesh const *pm) const {
  return flxup.GetSizeInBytes() + flxdn.GetSizeInBytes();
}

void Radiation::DumpRestartData(char *pdst) const {
  int offset = 0;

  std::memcpy(pdst + offset, flxup.data(), flxup.GetSizeInBytes());
  offset += flxup.GetSizeInBytes();
  std::memcpy(pdst + offset, flxdn.data(), flxdn.GetSizeInBytes());
  offset += flxdn.GetSizeInBytes();
}

//! \todo check me
size_t Radiation::LoadRestartData(char *psrc) {
  Application::Logger app("harp");
  app->Log("Radiation restarted");
  int offset = 0;

  std::memcpy(flxup.data(), psrc + offset, flxup.GetSizeInBytes());
  offset += flxup.GetSizeInBytes();
  std::memcpy(flxdn.data(), psrc + offset, flxdn.GetSizeInBytes());
  offset += flxdn.GetSizeInBytes();

  return offset;
}

namespace RadiationHelper {
bool real_close(Real num1, Real num2, Real tolerance) {
  return std::fabs(num1 - num2) <= tolerance;
}

std::pair<std::vector<Real>, std::vector<Real>> get_direction_grids(
    std::vector<Direction> const &dirs) {
  std::vector<Real> uphi;
  std::vector<Real> umu;

  for (auto &dir : dirs) {
    // find phi
    bool found = false;
    for (auto &phi : uphi)
      if (real_close(phi, dir.phi, 1.e-3)) {
        found = true;
        break;
      }
    if (!found) uphi.push_back(dir.phi);

    // find mu
    found = false;
    for (auto &mu : umu)
      if (real_close(mu, dir.mu, 1.e-3)) {
        found = true;
        break;
      }
    if (!found) umu.push_back(dir.mu);
  }

  std::sort(uphi.begin(), uphi.end());
  std::sort(umu.begin(), umu.end());

  return std::make_pair(uphi, umu);
}

Direction parse_radiation_direction(std::string_view str) {
  Direction ray;
  ray.phi = 0.;

  sscanf(str.data(), "(%lf,%lf)", &ray.mu, &ray.phi);
  ray.mu = cos(deg2rad(ray.mu));
  ray.phi = deg2rad(ray.phi);

  return ray;
}

std::vector<Direction> parse_radiation_directions(std::string str) {
  std::vector<std::string> dstr = Vectorize<std::string>(str.c_str());
  int nray = dstr.size();

  std::vector<Direction> ray(nray);

  auto jt = dstr.begin();
  for (auto it = ray.begin(); it != ray.end(); ++it, ++jt)
    *it = parse_radiation_direction(*jt);

  return ray;
}

uint64_t parse_radiation_flags(std::string str) {
  uint64_t flags = 0LL;
  std::stringstream msg;

  std::vector<std::string> dstr = Vectorize<std::string>(str.c_str(), " ,");
  for (int i = 0; i < dstr.size(); ++i) {
    if (dstr[i] == "time_dependent") {
      flags |= RadiationFlags::TimeDependent;
    } else if (dstr[i] == "broad_band") {
      flags |= RadiationFlags::BroadBand;
    } else if (dstr[i] == "stellar_beam") {
      flags |= RadiationFlags::StellarBeam;
    } else if (dstr[i] == "thermal_emission") {
      flags |= RadiationFlags::ThermalEmission;
    } else if (dstr[i] == "normalize") {
      flags |= RadiationFlags::Normalize;
    } else if (dstr[i] == "write_bin_radiance") {
      flags |= RadiationFlags::WriteBinRadiance;
    } else {
      msg << "flag: '" << dstr[i] << "' unrecognized" << std::endl;
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
