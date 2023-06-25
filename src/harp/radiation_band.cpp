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

// climath
#include <climath/core.h>

// application
#include <application/application.hpp>

// utils
#include <utils/fileio.hpp>
#include <utils/ndarrays.hpp>
#include <utils/vectorize.hpp>

// harp
#include "absorber.hpp"
#include "radiation.hpp"
#include "radiation_band.hpp"
#include "radiation_utils.hpp"  // readRadiationDirections

RadiationBand::RadiationBand(MeshBlock *pmb, ParameterInput *pin,
                             std::string name)
    : name_(name), bflags_(0LL), pmy_block_(pmb) {
  Application::Logger app("harp");
  app->Log("Initializing RadiationBand " + name_);
  std::stringstream msg;

  // parent radiation flags
  uint64_t rflags;
  set_radiation_flags(&rflags, pin->GetOrAddString("radiation", "flags", ""));

  // band flags
  if (pin->DoesParameterExist("radiation", name_ + ".flags")) {
    set_radiation_flags(&bflags_,
                        pin->GetString("radiation", name_ + ".flags"));
  }
  bflags_ |= rflags;

  // number of Legendre moments
  int npmom = pin->GetOrAddInteger("radiation", "npmom", 0);

  // name radiation band in the format of "min_wave max_wave nbins"
  std::string str = pin->GetString("radiation", name_);
  char default_file[80];
  snprintf(default_file, sizeof(default_file), "kcoeff.%s.nc", str.c_str());
  replaceChar(default_file, ' ', '-');

  std::vector<Real> val = Vectorize<Real>(str.c_str());
  if (val.size() != 3) {
    msg << "### FATAL ERROR in function RadiationBand::RadiationBand"
        << std::endl
        << "Length of '" << name_ << "' "
        << "must be 3.";
    ATHENA_ERROR(msg);
  }

  // set wavenumber and weights
  wmin_ = val[0];
  wmax_ = val[1];
  int num_bins = static_cast<int>(val[2]);
  if (num_bins < 1) {
    msg << "### FATAL ERROR in function RadiationBand::RadiationBand"
        << std::endl
        << "Length of some spectral band is not a positive number";
    ATHENA_ERROR(msg);
  }

  spec_.resize(num_bins);
  if (test(RadiationFlags::LineByLine)) {
    if (num_bins == 1) {
      if (wmin_ != wmax_) {
        msg << "### FATAL ERROR in function RadiationBand::RadiationBand"
            << std::endl
            << "The first spectrum must equal the last spectrum "
            << "if the length of the spectral band is 1.";
        ATHENA_ERROR(msg);
      }
      spec_[0].wav1 = spec_[0].wav2 = wmin_;
      spec_[0].wght = 1.;
    } else {
      Real dwave = (val[1] - val[0]) / (num_bins - 1);
      for (int i = 0; i < num_bins; ++i) {
        spec_[i].wav1 = spec_[i].wav2 = val[0] + dwave * i;
        spec_[i].wght = (i == 0) || (i == num_bins - 1) ? 0.5 * dwave : dwave;
      }
    }
  } else if (test(RadiationFlags::CorrelatedK)) {
    str = pin->GetString("radiation", name_ + ".gpoints");
    val = Vectorize<Real>(str.c_str(), ",");
    if (val.size() != num_bins) {
      msg << "### FATAL ERROR in function RadiationBand::RadiationBand"
          << std::endl
          << "Number of gpoints does not equal " << num_bins;
      ATHENA_ERROR(msg);
    }

    for (int i = 0; i < num_bins; ++i) spec_[i].wav1 = spec_[i].wav2 = val[i];

    str = pin->GetString("radiation", name_ + ".weights");
    val = Vectorize<Real>(str.c_str(), ",");
    if (val.size() != num_bins) {
      msg << "### FATAL ERROR in function RadiationBand::RadiationBand"
          << std::endl
          << "Number of weights does not equal " << num_bins;
      ATHENA_ERROR(msg);
    }

    for (int i = 0; i < num_bins; ++i) spec_[i].wght = val[i];
  } else {  // spectral bins
    Real dwave = (val[1] - val[0]) / num_bins;
    for (int i = 0; i < num_bins; ++i) {
      spec_[i].wav1 = val[0] + dwave * i;
      spec_[i].wav2 = val[0] + dwave * (i + 1);
      spec_[i].wght = 1.;
    }
  }

  // outgoing radiation direction (mu,phi) in degree
  if (pin->DoesParameterExist("radiation", name_ + ".outdir")) {
    str = pin->GetString("radiation", name_ + ".outdir");
    read_radiation_directions(&rayOutput_, str);
  } else if (pin->DoesParameterExist("radiation", "outdir")) {
    str = pin->GetString("radiation", "outdir");
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

  // absorbers
  str = pin->GetOrAddString("radiation", name + ".absorbers", "");
  std::vector<std::string> aname = Vectorize<std::string>(str.c_str());

  char astr[1024];
  for (int i = 0; i < aname.size(); ++i) {
    snprintf(astr, sizeof(astr), "%s.%s", name.c_str(), aname[i].c_str());
    std::string afile = pin->GetOrAddString("radiation", astr, default_file);
    addAbsorber(pin, name_, aname[i], afile);
  }

  // band parameters
  alpha_ = pin->GetOrAddReal("radiation", name + ".alpha", 0.);

  char buf[80];
  snprintf(buf, sizeof(buf), "%.2f - %.2f", wmin_, wmax_);
  app->Log("spectral range = " + std::string(buf));
  app->Log("number of spectral bins = " + std::to_string(num_bins));
}

RadiationBand::~RadiationBand() {
  for (size_t i = 0; i < absorbers.size(); ++i) delete absorbers[i];

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

// overide in the pgen file
void __attribute__((weak))
RadiationBand::addAbsorber(ParameterInput *pin, std::string bname,
                           std::string name, std::string file) {}

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
