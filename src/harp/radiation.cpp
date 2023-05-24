// C/C++ headers
#include <sstream>
#include <stdexcept>

// climath
#include <climath/core.h>

#include <athena/coordinates/coordinates.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// #include "../thermodynamics/thermodynamics.hpp"
// #include "../utils/utils.hpp"
#include <debugger/debugger.hpp>

#include "radiation.hpp"
#include "radiation_band.hpp"
#include "radiation_utils.hpp"  // setRadiationFlags

Real const Radiation::hPlanck = 6.63E-34;
Real const Radiation::hPlanck_cgs = 6.63E-27;
Real const Radiation::cLight = 3.E8;
Real const Radiation::cLight_cgs = 3.E10;
Real const Radiation::stefanBoltzmann = 5.670374419E-8;

Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin)
    : rflags_(0LL), pcoord_(pmb->pcoord) {
  pdebug->Enter("Radiation");

  // radiation flags
  set_radiation_flags(&rflags_, pin->GetOrAddString("radiation", "flags", ""));

  // distance to parent star
  stellarDistance_au_ = pin->GetOrAddReal("radiation", "distance_au", 1.);
  pdebug->Message("stellar distance", stellarDistance_au_);

  // radiation bands
  int bid = 1;
  readRadiationBands(pmb, pin, bid);

  // incoming radiation direction (mu,phi) in degree
  std::string str = pin->GetOrAddString("radiation", "indir", "(0.,0.)");
  read_radiation_directions(rayInput_, str);

  // output radiance
  int nout = getNumOutgoingRays();
  if (nout > 0) {
    radiance.NewAthenaArray(nout, pmb->ncells3, pmb->ncells2);
  }

  // time control
  cooldown_ = pin->GetOrAddReal("radiation", "dt", 0.);
  current_ = 0.;

  planet_ = new CelestrialBody(pin);
  pdebug->Leave();
}

Radiation::~Radiation() {
  for (size_t i = 0; i < bands.size(); ++i) {
    delete bands[i];
  }
  delete planet_;
}

void Radiation::calculateRadiativeFlux(AthenaArray<Real> &flxup,
                                       AthenaArray<Real> &flxdn, Real time,
                                       int k, int j, int il, int iu) {
  pdebug->Call("Radiation::calculateRadiativeFlux");
  Real dist = stellarDistance_au_;

  int idx = 0;
  for (auto p : bands) {
    Direction ray;
    if (p->test(RadiationFlags::Dynamic)) {
      planet_->ParentZenithAngle(&ray.mu, &ray.phi, time, pcoord_->x2v(j),
                                 pcoord_->x3v(k));
      dist = planet_->ParentDistanceInAu(time);
    } else {
      ray = rayInput_[0];
    }

    // iu ~= ie + 1
    AthenaArray<Real> bflxup, bflxdn;
    bflxup.InitWithShallowSlice(flxup, 3, idx, p->getNumOutgoingRays());
    bflxdn.InitWithShallowSlice(flxdn, 3, idx, p->getNumOutgoingRays());

    p->setSpectralProperties(k, j, il - NGHOST, iu + NGHOST - 1);
    p->calculateBandFlux(bflxup, bflxdn, ray, dist, k, j, il, iu);
    idx++;
  }

  pdebug->Leave();
}

void Radiation::calculateRadiance(AthenaArray<Real> &radiance, Real time, int k,
                                  int j, int il, int iu) {
  pdebug->Call("Radiation::calculateRadiance");
  Real dist = stellarDistance_au_;

  Direction ray;
  if (test(RadiationFlags::Dynamic)) {
    planet_->ParentZenithAngle(&ray.mu, &ray.phi, time, pcoord_->x2v(j),
                               pcoord_->x3v(k));
    dist = planet_->ParentDistanceInAu(time);
  }

  int idx = 0;
  for (auto p : bands) {
    if (p->getNumOutgoingRays() == 0) continue;

    // iu ~= ie + 1
    AthenaArray<Real> brad;
    brad.InitWithShallowSlice(radiance, 3, idx, p->getNumOutgoingRays());

    p->setSpectralProperties(k, j, il - NGHOST, iu + NGHOST - 1);
    p->calculateBandRadiance(brad, ray, dist, k, j, il, iu);
    idx += p->getNumOutgoingRays();
  }

  pdebug->Leave();
}

void Radiation::addRadiativeFlux(Hydro *phydro, int k, int j, int il,
                                 int iu) const {
  pdebug->Call("Radiation::addRadiativeFlux");

  // x1-flux divergence
  for (size_t b = 0; b < bands.size(); ++b) {
#pragma omp simd
    for (int i = il; i <= iu; ++i)
      phydro->flux[X1DIR](IEN, k, j, i) +=
          flxup(b, k, j, i) - flxdn(b, k, j, i);
  }

  pdebug->Leave();
}

void Radiation::readRadiationBands(MeshBlock *pmb, ParameterInput *pin,
                                   int &bid) {
  char name[80];
  while (true) {
    snprintf(name, 80, "b%d", bid);
    if (!pin->DoesParameterExist("radiation", name)) break;
    RadiationBand *p = new RadiationBand(pmb, pin, name);
    bands.push_back(p);
    bid++;
  }

  if (pin->DoesParameterExist("radiation", "bandsfile")) {
    ParameterInput *pin_next = new ParameterInput;
    IOWrapper infile;
    infile.Open(pin->GetString("radiation", "bandsfile").c_str(),
                IOWrapper::FileMode::read);
    pin_next->LoadFromFile(infile);
    infile.Close();
    InputBlock *pblock = pin->GetPtrToBlock("radiation");
    InputLine *pline = pblock->pline;

    // remove the bandsfile line
    while (pline->pnext != nullptr) {
      if (pline->pnext->param_name == "bandsfile") {
        InputLine *pnext = pline->pnext->pnext;
        delete pline->pnext;
        pline->pnext = pnext;
        continue;
      }
      pline = pline->pnext;
    }

    // get the first line of the current input in block radiation
    pline = pin_next->GetPtrToBlock("radiation")->pline;

    // copy the current lines into the main input
    while (pline != nullptr) {
      pin->AddParameter(pblock, pline->param_name, pline->param_value,
                        pline->param_comment);
      pline = pline->pnext;
    }
    readRadiationBands(pmb, pin, bid);
    delete pin_next;
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
