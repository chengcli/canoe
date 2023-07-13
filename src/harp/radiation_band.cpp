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
#include "rt_solvers.hpp"

RadiationBand::RadiationBand(MeshBlock *pmb, ParameterInput *pin,
                             YAML::Node &node, std::string myname)
    : name_(myname), bflags_(0LL), pmy_block_(pmb) {
  Application::Logger app("harp");
  app->Log("Initialize RadiationBand " + name_);

  if (!node[name_]) {
    throw NotFoundError("RadiationBand", name_);
  }

  auto my = node[name_];

  if (my["parameters"]) params_ = ToParameterMap(my["parameters"]);

  type_ = my["type"] ? my["type"].as<std::string>() : "unknown";

  // number of Legendre moments
  int npmom = my["moments"] ? my["moments"].as<int>() : 1;

  // set wavenumber and weights
  if (my["wavenumber-range"]) {
    setWavenumberRange(my);
  } else if (my["wavenumbers"]) {
    setWavenumberGrid(my);
  } else if (my["wavelength-range"]) {
    setWavelengthRange(my);
  } else if (my["wavelengths"]) {
    setWavelengthGrid(my);
  } else if (my["frequency-Range"]) {
    setFrequencyRange(my);
  } else if (my["frequencies"]) {
    setFrequencyGrid(my);
  } else {
    throw NotFoundError("RadiationBand", "Spectral range");
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

  tau_.NewAthenaArray(spec_.size(), ncells1);
  tau_.ZeroClear();

  ssa_.NewAthenaArray(spec_.size(), ncells1);
  ssa_.ZeroClear();

  toa_.NewAthenaArray(spec_.size(), rayOutput_.size());
  toa_.ZeroClear();

  pmom_.NewAthenaArray(spec_.size(), ncells1, npmom + 1);
  pmom_.ZeroClear();

  flxup_.NewAthenaArray(spec_.size(), ncells1);
  flxup_.ZeroClear();

  flxdn_.NewAthenaArray(spec_.size(), ncells1);
  flxdn_.ZeroClear();

  // band properties
  btau.NewAthenaArray(ncells3, ncells2, ncells1);
  bssa.NewAthenaArray(ncells3, ncells2, ncells1);
  bpmom.NewAthenaArray(npmom + 1, ncells3, ncells2, ncells1);

  //! \note btoa, bflxup, bflxdn are shallow slices to Radiation variables

  // add absorbers
  if (my["opacity"]) {
    for (auto aname : my["opacity"]) {
      bool found = false;
      for (auto absorber : node["opacity-sources"]) {
        if (type_ + "-" + aname.as<std::string>() ==
            absorber["name"].as<std::string>()) {
          AddAbsorber(pin, absorber);
          found = true;
          break;
        }
      }

      if (!found) {
        throw NotFoundError("RadiationBand",
                            "Opacity " + aname.as<std::string>());
      }
    }
  } else {
    throw NotFoundError("RadiationBand", "Band " + name_ + " opacity");
  }

  // set rt solver
  if (my["rt-solver"]) {
    if (my["rt-solver"].as<std::string>() == "Lambert") {
      psolver = std::make_shared<RTSolverLambert>(this);
#ifdef RT_DISORT
    } else if (my["rt-solver"].as<std::string>() == "Disort") {
      psolver = std::make_shared<RTSolverDisort>(this);
#endif
    } else {
      throw RuntimeError("RadiationBand", my["rt-solver"].as<std::string>());
    }
  } else {
    psolver = std::make_shared<RTSolver>(this, "Null");
  }

  char buf[80];
  snprintf(buf, sizeof(buf), "%.2f - %.2f", wmin_, wmax_);
  app->Log("Spectral range", buf);
  app->Log("Number of spectral bins", spec_.size());
}

RadiationBand::~RadiationBand() {
  Application::Logger app("harp");
  app->Log("Destroy RadiationBand " + name_);
}

void RadiationBand::setWavenumberRange(YAML::Node &my) {
  wmin_ = my["wavenumber-range"][0].as<Real>();
  wmax_ = my["wavenumber-range"][1].as<Real>();
  if (wmin_ > wmax_) {
    throw RuntimeError("setWavenumberRange", "wmin > wmax");
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
        "RadiationBand",
        "either 'resolution' or 'num-bins' or 'gpoints' must be defined");
  }
}

void RadiationBand::setWavenumberGrid(YAML::Node &my) {
  throw NotImplementedError("setWavenumberGrid");
}

void RadiationBand::setFrequencyRange(YAML::Node &my) {
  throw NotImplementedError("setFrequencyRange");
}

void RadiationBand::setFrequencyGrid(YAML::Node &my) {
  std::vector<Real> freqs = my["frequencies"].as<std::vector<Real>>();

  wmin_ = freqs.front();
  wmax_ = freqs.back();

  if (wmin_ > wmax_) {
    throw RuntimeError("setWavenumberRange", "wmin > wmax");
  }

  spec_.resize(freqs.size());

  for (int i = 0; i < spec_.size(); ++i) {
    spec_[i].wav1 = freqs[i];
    spec_[i].wav2 = freqs[i];
    spec_[i].wght = 1.;
  }
}

void RadiationBand::setWavelengthRange(YAML::Node &my) {
  throw NotImplementedError("setWavelengthRange");
}

void RadiationBand::setWavelengthGrid(YAML::Node &my) {
  throw NotImplementedError("setWavelengthGrid");
}

Absorber *RadiationBand::GetAbsorber(std::string const &name) {
  for (auto &absorber : absorbers_) {
    if (absorber->GetName() == name) {
      return absorber.get();
    }
  }

  throw NotFoundError("RadiationBand", "Absorber " + name);

  return nullptr;
}

void RadiationBand::WriteBinRadiance(OutputParameters const *pout) const {
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

void RadiationBand::AddAbsorber(ParameterInput *pin, YAML::Node &node) {
  if (strcmp(PLANET, "Jupiter") == 0 || strcmp(PLANET, "Saturn") == 0) {
    addAbsorberGiants(pin, node);
  } else if (strcmp(PLANET, "Uranus") == 0 || strcmp(PLANET, "Neptune") == 0) {
    addAbsorberGiants(pin, node);
  } else if (strcmp(PLANET, "Earth") == 0) {
    addAbsorberEarth(pin, node);
  } else if (strcmp(PLANET, "Mars") == 0) {
    addAbsorberMars(pin, node);
  } else if (strcmp(PLANET, "Venus") == 0) {
    addAbsorberVenus(pin, node);
  }
}
