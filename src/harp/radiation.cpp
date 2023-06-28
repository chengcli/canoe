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

Real const Radiation::hPlanck = 6.63E-34;
Real const Radiation::hPlanck_cgs = 6.63E-27;
Real const Radiation::cLight = 3.E8;
Real const Radiation::cLight_cgs = 3.E10;
Real const Radiation::stefanBoltzmann = 5.670374419E-8;

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
  if (pin->DoesParameterExist("radiation", "bandsfile"))
    PopulateRadiationBands(pin);

  // incoming radiation direction (mu,phi) in degree
  std::string str = pin->GetOrAddString("radiation", "indir", "(0.,0.)");
  read_radiation_directions(&rayInput_, str);

  // output radiance
  int nout = getNumOutgoingRays();
  if (nout > 0) {
    radiance.NewAthenaArray(nout, pmb->ncells3, pmb->ncells2);
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

RadiationBand *Radiation::GetBand(std::string const &name) {
  for (auto &band : bands) {
    if (band->GetName() == name) {
      return band.get();
    }
  }

  throw NotFoundError("GetBand", "Band " + name);
}

void Radiation::PopulateRadiationBands(ParameterInput *pin) {
  Application::Logger app("harp");
  app->Log("Populate Radiation bands");

  std::string filename = pin->GetString("radiation", "bandsfile");
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
    bands.push_back(std::move(p));
  }
}

void Radiation::calculateRadiativeFlux(AthenaArray<Real> *flxup,
                                       AthenaArray<Real> *flxdn, Real time,
                                       int k, int j, int il, int iu) {
  Application::Logger app("harp");
  app->Log("CalculateRadiativeFlux");
  Real dist = stellarDistance_au_;

  int idx = 0;
  Coordinates *pcoord = pmy_block_->pcoord;
  for (auto p : bands) {
    Direction ray;
    if (p->test(RadiationFlags::Dynamic)) {
      planet_->ParentZenithAngle(&ray.mu, &ray.phi, time, pcoord->x2v(j),
                                 pcoord->x3v(k));
      dist = planet_->ParentDistanceInAu(time);
    } else {
      ray = rayInput_[0];
    }

    // iu ~= ie + 1
    AthenaArray<Real> bflxup, bflxdn;
    bflxup.InitWithShallowSlice(*flxup, 3, idx, p->getNumOutgoingRays());
    bflxdn.InitWithShallowSlice(*flxdn, 3, idx, p->getNumOutgoingRays());

    p->SetSpectralProperties(k, j, il - NGHOST, iu + NGHOST - 1);
    p->calculateBandFlux(&bflxup, &bflxdn, ray, dist, k, j, il, iu);
    idx++;
  }
}

void Radiation::calculateRadiance(AthenaArray<Real> *radiance, Real time, int k,
                                  int j, int il, int iu) {
  Application::Logger app("harp");
  app->Log("CalculateRadiance");
  Real dist = stellarDistance_au_;
  Coordinates *pcoord = pmy_block_->pcoord;

  Direction ray;
  if (test(RadiationFlags::Dynamic)) {
    planet_->ParentZenithAngle(&ray.mu, &ray.phi, time, pcoord->x2v(j),
                               pcoord->x3v(k));
    dist = planet_->ParentDistanceInAu(time);
  }

  int idx = 0;
  for (auto p : bands) {
    if (p->getNumOutgoingRays() == 0) continue;

    // iu ~= ie + 1
    AthenaArray<Real> brad;
    brad.InitWithShallowSlice(*radiance, 3, idx, p->getNumOutgoingRays());

    p->SetSpectralProperties(k, j, il - NGHOST, iu + NGHOST - 1);
    p->calculateBandRadiance(&brad, ray, dist, k, j, il, iu);
    idx += p->getNumOutgoingRays();
  }
}

void Radiation::addRadiativeFlux(Hydro *phydro, int k, int j, int il,
                                 int iu) const {
  Application::Logger app("harp");
  app->Log("AddRadiativeFlux");

  // x1-flux divergence
  for (size_t b = 0; b < bands.size(); ++b) {
#pragma omp simd
    for (int i = il; i <= iu; ++i)
      phydro->flux[X1DIR](IEN, k, j, i) +=
          flxup(b, k, j, i) - flxdn(b, k, j, i);
  }
}

size_t Radiation::getNumOutgoingRays() const {
  size_t num = 0;
  for (auto p : bands) {
    num += p->getNumOutgoingRays();
  }
  return num;
}

size_t Radiation::getRestartDataSizeInBytes() const {
  return flxup.GetSizeInBytes() + flxdn.GetSizeInBytes();
}

size_t Radiation::dumpRestartData(char *pdst) const {
  int offset = 0;

  std::memcpy(pdst + offset, flxup.data(), flxup.GetSizeInBytes());
  offset += flxup.GetSizeInBytes();
  std::memcpy(pdst + offset, flxdn.data(), flxdn.GetSizeInBytes());
  offset += flxdn.GetSizeInBytes();

  return getRestartDataSizeInBytes();
}

size_t Radiation::loadRestartData(char *psrc) {
  int offset = 0;

  std::memcpy(flxup.data(), psrc + offset, flxup.GetSizeInBytes());
  offset += flxup.GetSizeInBytes();
  std::memcpy(flxdn.data(), psrc + offset, flxdn.GetSizeInBytes());
  offset += flxdn.GetSizeInBytes();

  return getRestartDataSizeInBytes();
}
