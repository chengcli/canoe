// C/C++ headers
#include <fstream>
#include <sstream>
#include <stdexcept>

// external
#include <yaml-cpp/yaml.h>

// climath
#include <climath/core.h>

// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// harp
#include "radiation.hpp"
#include "radiation_band.hpp"
#include "radiation_utils.hpp"  // setRadiationFlags
#include "rt_solvers.hpp"

Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin)
    : rflags_(0LL), pmy_block_(pmb) {
  Application::Logger app("harp");
  app->Log("Initialize Radiation");

  // radiation flags
  set_radiation_flags(&rflags_, pin->GetOrAddString("radiation", "flags", ""));

  // distance to parent star
  stellarDistance_au_ = pin->GetOrAddReal("radiation", "distance_au", 1.);
  app->Log("Stellar distance = " + std::to_string(stellarDistance_au_) + " au");

  // radiation bands
  if (pin->DoesParameterExist("radiation", "control_file"))
    PopulateRadiationBands(pin);

  // incoming radiation direction (mu,phi) in degree
  std::string str = pin->GetOrAddString("radiation", "indir", "(0.,0.)");
  read_radiation_directions(&rayInput_, str);

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
  cooldown_ = pin->GetOrAddReal("radiation", "dt", 0.);
  current_ = 0.;

  planet_ = std::make_unique<CelestrialBody>(pin);
}

Radiation::~Radiation() {
  Application::Logger app("harp");
  app->Log("Destroy Radiation");
}

RadiationBand *Radiation::GetBand(std::string const &name) const {
  for (auto &band : bands_) {
    if (band->GetName() == name) {
      return band.get();
    }
  }

  throw NotFoundError("GetBand", "Band " + name);
}

void Radiation::PopulateRadiationBands(ParameterInput *pin) {
  Application::Logger app("harp");
  app->Log("Populate Radiation bands");

  std::string filename = pin->GetString("radiation", "control_file");

  std::ifstream stream(filename);
  if (stream.good() == false) {
    app->Error("Cannot open radiation bands file: " + filename);
  }
  YAML::Node node = YAML::Load(stream);

  if (!node["opacity-sources"]) {
    throw NotFoundError("PopulateRadiationBand", "opacity-sources");
  }

  for (auto bname : node["bands"]) {
    auto p = std::make_unique<RadiationBand>(pmy_block_, pin, node,
                                             bname.as<std::string>());
    bands_.push_back(std::move(p));
  }
}

void Radiation::CalRadiativeFlux(Real time, int k, int j, int il, int iu) {
  Application::Logger app("harp");
  app->Log("CalculateRadiativeFlux");
  Real dist = stellarDistance_au_;

  Coordinates *pcoord = pmy_block_->pcoord;
  Direction ray;
  if (test(RadiationFlags::Dynamic)) {
    planet_->ParentZenithAngle(&ray.mu, &ray.phi, time, pcoord->x2v(j),
                               pcoord->x3v(k));
    dist = planet_->ParentDistanceInAu(time);
  } else {
    ray = rayInput_[0];
  }

  for (auto &p : bands_) {
    // iu ~= ie + 1
    p->SetSpectralProperties(k, j, il - NGHOST, iu + NGHOST - 1);
    p->psolver->CalBandFlux(ray, dist, k, j, il, iu);
  }
}

void Radiation::CalRadiance(Real time, int k, int j, int il, int iu) {
  Application::Logger app("harp");
  app->Log("CalculateRadiance");
  Real dist = stellarDistance_au_;

  Coordinates *pcoord = pmy_block_->pcoord;
  Direction ray;
  if (test(RadiationFlags::Dynamic)) {
    planet_->ParentZenithAngle(&ray.mu, &ray.phi, time, pcoord->x2v(j),
                               pcoord->x3v(k));
    dist = planet_->ParentDistanceInAu(time);
  } else {
    ray = rayInput_[0];
  }

  for (auto &p : bands_) {
    // iu ~= ie + 1
    p->SetSpectralProperties(k, j, il - NGHOST, iu + NGHOST - 1);
    p->psolver->CalBandRadiance(ray, dist, k, j, il, iu);
  }
}

void Radiation::AddRadiativeFlux(Hydro *phydro, int k, int j, int il,
                                 int iu) const {
  Application::Logger app("harp");
  app->Log("AddRadiativeFlux");

  // x1-flux divergence
  for (size_t b = 0; b < bands_.size(); ++b) {
#pragma omp simd
    for (int i = il; i <= iu; ++i)
      phydro->flux[X1DIR](IEN, k, j, i) +=
          flxup(b, k, j, i) - flxdn(b, k, j, i);
  }
}

size_t Radiation::GetNumOutgoingRays() const {
  size_t num = 0;
  for (auto &p : bands_) {
    num += p->GetNumOutgoingRays();
  }
  return num;
}

size_t Radiation::GetRestartDataSizeInBytes() const {
  return flxup.GetSizeInBytes() + flxdn.GetSizeInBytes();
}

size_t Radiation::DumpRestartData(char *pdst) const {
  int offset = 0;

  std::memcpy(pdst + offset, flxup.data(), flxup.GetSizeInBytes());
  offset += flxup.GetSizeInBytes();
  std::memcpy(pdst + offset, flxdn.data(), flxdn.GetSizeInBytes());
  offset += flxdn.GetSizeInBytes();

  return GetRestartDataSizeInBytes();
}

size_t Radiation::LoadRestartData(char *psrc) {
  int offset = 0;

  std::memcpy(flxup.data(), psrc + offset, flxup.GetSizeInBytes());
  offset += flxup.GetSizeInBytes();
  std::memcpy(flxdn.data(), psrc + offset, flxdn.GetSizeInBytes());
  offset += flxdn.GetSizeInBytes();

  return GetRestartDataSizeInBytes();
}
