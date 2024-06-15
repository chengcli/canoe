// cantera
#include <cantera/base/yaml.h>

// athena
#include <athena/athena.hpp>

// application
#include <application/exceptions.hpp>

// harp
#include "radiation.hpp"
#include "spectral_grid.hpp"

SpectralGridBase::SpectralGridBase(YAML::Node const& my)
    : unit_(RadiationHelper::parse_unit_with_default(my)) {}

RegularSpacingSpectralGrid::RegularSpacingSpectralGrid(YAML::Node const& my)
    : SpectralGridBase(my) {
  auto wpair = RadiationHelper::parse_wave_range(my);

  Real wmin = wpair.first;
  Real wmax = wpair.second;
  int num_bins;

  if (wmin == wmax) {
    num_bins = 1;
    spec.resize(num_bins);
    spec[0].wav1 = spec[0].wav2 = wmin;
    spec[0].wght = 1.;
  } else if (my["resolution"]) {
    Real dwave = my["resolution"].as<Real>();
    num_bins = static_cast<int>((wmax - wmin) / dwave) + 1;
    spec.resize(num_bins);
    for (int i = 0; i < num_bins; ++i) {
      spec[i].wav1 = spec[i].wav2 = wmin + dwave * i;
      spec[i].wght = (i == 0) || (i == num_bins - 1) ? 0.5 * dwave : dwave;
    }
  } else if (my["num-bins"]) {
    num_bins = my["num-bins"].as<int>();
    Real dwave = static_cast<Real>(1. * (wmax - wmin) / num_bins);
    spec.resize(num_bins);
    for (int i = 0; i < num_bins; ++i) {
      spec[i].wav1 = wmin + dwave * i;
      spec[i].wav2 = spec[i].wav1 + dwave;
      spec[i].wght = 1.;
    }
  } else {
    throw NotFoundError("RegularSpacingSpectralGrid",
                        "either 'resolution' or 'num-bins' must be defined");
  }

  // normalize weights to add up to 1
  Real total_wght = 0;

  for (int i = 0; i < num_bins; ++i) total_wght += spec[i].wght;

  for (int i = 0; i < num_bins; ++i) spec[i].wght /= total_wght;
}

CustomSpacingSpectralGrid::CustomSpacingSpectralGrid(YAML::Node const& my)
    : SpectralGridBase(my) {
  char str[80];
  snprintf(str, sizeof(str), "%s-points", unit_.c_str());

  if (!my[str]) {
    throw NotFoundError("CustomSpacingSpectralGrid", str);
  }

  std::vector<Real> wavs = my[str].as<std::vector<Real>>();
  std::sort(wavs.begin(), wavs.end());

  Real wmin = wavs.front();
  Real wmax = wavs.back();
  int num_bins = wavs.size();

  if (wmin > wmax) {
    throw RuntimeError("CustomSpacingSpectralGrid", "wmin > wmax");
  }

  spec.resize(num_bins);
  for (int i = 0; i < num_bins; ++i) {
    spec[i].wav1 = wavs[i];
    spec[i].wav2 = wavs[i];
    spec[i].wght = 1.;
  }
}

// delayed setting weights
CKTableSpectralGrid::CKTableSpectralGrid(YAML::Node const& my)
    : SpectralGridBase(my) {}
