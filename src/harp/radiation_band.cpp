// C/C++
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

// athena
#include <athena/mesh/mesh.hpp>
#include <athena/outputs/outputs.hpp>
#include <athena/parameter_input.hpp>
#include <athena/utils/utils.hpp>

// canoe
#include <configure.hpp>

// climath
#include <climath/core.h>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// utils
#include <utils/fileio.hpp>
#include <utils/ndarrays.hpp>
#include <utils/parameter_map.hpp>
#include <utils/vectorize.hpp>

// harp
#include "absorber.hpp"
#include "radiation.hpp"
#include "radiation_band.hpp"
#include "radiation_utils.hpp"  // readRadiationDirections

RadiationBand::RadiationBand(MeshBlock *pmb, ParameterInput *pin,
                             YAML::Node &node, std::string myname)
    : name_(myname),
      bflags_(0LL),
      pmy_block_(pmb),
      params_(ToParameterMap(node["parameters"])) {
  Application::Logger app("harp");
  app->Log("Initialize RadiationBand " + name_);

  if (!node[name_]) {
    throw NotFoundError("RadiationBand " + name_);
  }

  auto my = node[name_];

  // number of Legendre moments
  int npmom = my["moments"] ? my["moments"].as<int>() : 1;

  // set wavenumber and weights
  wmin_ = my["wavenumber-range"][0].as<Real>();
  wmax_ = my["wavenumber-range"][1].as<Real>();
  if (wmin_ >= wmax_) {
    app->Error("Wavenumber range is not valid " + std::to_string(wmin_) + " " +
               std::to_string(wmax_));
  }

  int num_bins = 1;
  if (wmin_ == wmax_) {
    num_bins = 1;
    spec_.resize(num_bins);
    spec_[0].wav1 = spec_[0].wav2 = wmin_;
    spec_[0].wght = 1.;
  } else if (my["resolution"]) {
    Real dwave = my["resolution"].as<Real>();
    num_bins = static_cast<int>((wmax_ - wmin_) / dwave) + 1;
    spec_.resize(num_bins);
    for (int i = 0; i < num_bins; ++i) {
      spec_[i].wav1 = spec_[i].wav2 = wmin_ + dwave * i;
      spec_[i].wght = (i == 0) || (i == num_bins - 1) ? 0.5 * dwave : dwave;
    }
  } else if (my["num-bins"]) {
    Real dwave = static_cast<Real>(1. * (wmax_ - wmin_) / num_bins);
    num_bins = my["num-bins"].as<int>();
    spec_.resize(num_bins);
    for (int i = 0; i < num_bins; ++i) {
      spec_[i].wav1 = wmin_ + dwave * i;
      spec_[i].wav2 = spec_[i].wav1 + dwave;
      spec_[i].wght = 1.;
    }
  } else if (my["gpoints"]) {
    num_bins = my["gpoints"].size();
    spec_.resize(num_bins);
    for (int i = 0; i < num_bins; ++i) {
      spec_[i].wav1 = wmin_;
      spec_[i].wav2 = wmax_;
      spec_[i].wght = my["gpoints"][i].as<Real>();
    }
  } else {
    throw NotFoundError(
        "either 'resolution' or 'num-bins' or 'gpoints' must be defined");
  }

  // outgoing radiation direction (mu,phi) in degree
  if (pin->DoesParameterExist("radiation", name_ + ".outdir")) {
    auto str = pin->GetString("radiation", name_ + ".outdir");
    read_radiation_directions(&rayOutput_, str);
  } else if (pin->DoesParameterExist("radiation", "outdir")) {
    auto str = pin->GetString("radiation", "outdir");
    read_radiation_directions(&rayOutput_, str);
  }

  // allocate memory
  int ncells1 = pmb->ncells1;
  int ncells2 = pmb->ncells2;
  int ncells3 = pmb->ncells3;

  // spectral properties
  tem_.NewAthenaArray(ncells1);
  temf_.NewAthenaArray(ncells1 + 1);

  tau_.NewAthenaArray(num_bins, ncells1);
  tau_.ZeroClear();

  ssa_.NewAthenaArray(num_bins, ncells1);
  ssa_.ZeroClear();

  toa_.NewAthenaArray(num_bins, rayOutput_.size());
  toa_.ZeroClear();

  pmom_.NewAthenaArray(num_bins, ncells1, npmom + 1);
  pmom_.ZeroClear();

  // band properties
  btau.NewAthenaArray(ncells3, ncells2, ncells1);
  bssa.NewAthenaArray(ncells3, ncells2, ncells1);
  bpmom.NewAthenaArray(npmom + 1, ncells3, ncells2, ncells1);

  // add absorbers
  if (my["opacity"]) {
    for (auto aname : my["opacity"]) {
      bool found = false;
      for (auto absorber : node["opacity-sources"]) {
        if (aname.as<std::string>() == absorber["name"].as<std::string>()) {
          AddAbsorber(pin, name_, absorber);
          found = true;
          break;
        }
      }

      if (!found) {
        throw NotFoundError("Opacity " + aname.as<std::string>());
      }
    }
  } else {
    throw NotFoundError("Band " + name_ + " opacity");
  }

  char buf[80];
  snprintf(buf, sizeof(buf), "%.2f - %.2f", wmin_, wmax_);
  app->Log("Spectral range = " + std::string(buf));
  app->Log("Number of spectral bins = " + std::to_string(num_bins));
}

RadiationBand::~RadiationBand() {
  Application::Logger app("harp");
  app->Log("Destroy RadiationBand " + name_);

#ifdef RT_DISORT
  free_disort();
#endif
}

void RadiationBand::writeBinRadiance(OutputParameters const *pout) const {
  if (!test(RadiationFlags::WriteBinRadiance)) return;

  char fname[80], number[6];
  snprintf(number, sizeof(number), "%05d", pout->file_number);
  snprintf(fname, sizeof(fname), "%s.radiance.%s.txt", name_.c_str(), number);
  FILE *pfile = fopen(fname, "w");

  fprintf(pfile, "# Bin Radiances of Band %s: %.3g - %.3g\n", name_.c_str(),
          wmin_, wmax_);
  fprintf(pfile, "# Ray output size: %lu\n", rayOutput_.size());

  fprintf(pfile, "# Polar angles: ");
  for (size_t j = 0; j < rayOutput_.size(); ++j) {
    fprintf(pfile, "%.3f", rad2deg(acos(rayOutput_[j].mu)));
  }
  fprintf(pfile, "\n");

  fprintf(pfile, "# Azimuthal angles: ");
  for (size_t j = 0; j < rayOutput_.size(); ++j) {
    fprintf(pfile, "%.3f", rad2deg(rayOutput_[j].phi));
  }
  fprintf(pfile, "\n");

  fprintf(pfile, "#%12s%12s", "wave1", "wave2");
  for (size_t j = 0; j < rayOutput_.size(); ++j) {
    fprintf(pfile, "%12s%lu", "Radiance", j + 1);
  }
  fprintf(pfile, "%12s\n", "weight");

  for (size_t i = 0; i < spec_.size(); ++i) {
    fprintf(pfile, "%13.3g%12.3g", spec_[i].wav1, spec_[i].wav2);
    for (size_t j = 0; j < rayOutput_.size(); ++j) {
      fprintf(pfile, "%12.3f", toa_(i, j));
    }
    if (test(RadiationFlags::Normalize) && (wmax_ != wmin_)) {
      fprintf(pfile, "%12.3g\n", spec_[i].wght / (wmax_ - wmin_));
    } else {
      fprintf(pfile, "%12.3g\n", spec_[i].wght);
    }
  }

  fclose(pfile);
}

void RadiationBand::AddAbsorber(ParameterInput *pin, std::string bname,
                                YAML::Node &node) {
  if (PLANET == "Jupiter" || PLANET == "Saturn") {
    addAbsorberGiants(pin, bname, node);
  } else if (PLANET == "Uranus" || PLANET == "Neptune") {
    addAbsorberGiants(pin, bname, node);
  } else if (PLANET == "Earth") {
    addAbsorberEarth(pin, bname, node);
  } else if (PLANET == "Mars") {
    addAbsorberMars(pin, bname, node);
  } else if (PLANET == "Venus") {
    addAbsorberVenus(pin, bname, node);
  } else {
    throw NotFoundError(PLANET);
  }
}

// overide in rtsolver folder
void __attribute__((weak))
RadiationBand::calculateBandFlux(AthenaArray<Real> *flxup,
                                 AthenaArray<Real> *flxdn,
                                 Direction const &rayInput, Real dist, int k,
                                 int j, int il, int iu) {}

// overide in rtsolver folder
void __attribute__((weak))
RadiationBand::calculateBandRadiance(AthenaArray<Real> *radiance,
                                     Direction const &rayInput, Real dist,
                                     int k, int j, int il, int iu) {}
