// external
#include <yaml-cpp/yaml.h>

// athena
#include <athena/athena.hpp>

// application
#include <application/exceptions.hpp>

// harp
#include "spectral_grid.hpp"

std::pair<Real, Real> SpectralGridBase::ReadRangeFrom(YAML::Node const& my) {
  char str[80];
  snprintf(str, sizeof(str), "%s-%s", unit_type.c_str(), "range");

  if (!my[str]) {
    throw NotFoundError("SpectralGridBase", str);
  }

  /// wavenumber-range, wavelength-range, frequency-range, etc
  Real wmin = my[str][0].as<Real>();
  Real wmax = my[str][1].as<Real>();
  if (wmin > wmax) {
    throw RuntimeError("SpectralGridBase", "wmin > wmax");
  }

  return std::make_pair(wmin, wmax);
}

RegularSpacingSpectralGrid::RegularSpacingSpectralGrid(YAML::Node const& my) {
  std::string units = my["units"] ? my["units"].as<std::string>() : "cm-1";

  if (units == "cm-1") {
    unit_type = "wavenumber";
  } else if (units == "um") {
    unit_type = "wavelength";
  } else if (units == "GHz") {
    unit_type = "frequency";
  } else {
    throw RuntimeError("RegularSpacingSpectralGrid",
                       "unknown spectral unit type");
  }

  auto&& wpair = ReadRangeFrom(my);
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
    Real dwave = static_cast<Real>(1. * (wmax - wmin) / num_bins);
    num_bins = my["num-bins"].as<int>();
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
}

CustomSpacingSpectralGrid::CustomSpacingSpectralGrid(YAML::Node const& my) {
  char str[80];
  snprintf(str, sizeof(str), "%s-points", unit_type.c_str());

  if (!my[str]) {
    throw NotFoundError("SpectralGridBase", str);
  }

  std::vector<Real> wavs = my[str].as<std::vector<Real>>();
  std::sort(wavs.begin(), wavs.end());

  Real wmin = wavs.front();
  Real wmax = wavs.back();
  int num_bins = wavs.size();

  if (wmin > wmax) {
    throw RuntimeError("SpectralGridBase", "wmin > wmax");
  }

  spec.resize(num_bins);
  for (int i = 0; i < num_bins; ++i) {
    spec[i].wav1 = wavs[i];
    spec[i].wav2 = wavs[i];
    spec[i].wght = 1.;
  }
}

CorrelatedKTableSpectralGrid::CorrelatedKTableSpectralGrid(
    YAML::Node const& my) {
  auto&& wpair = ReadRangeFrom(my);
  Real wmin = wpair.first;
  Real wmax = wpair.second;

  if (my["g-points"]) {
    int num_bins = my["g-points"].size();
    spec.resize(num_bins);
    for (int i = 0; i < num_bins; ++i) {
      spec[i].wav1 = wmin;
      spec[i].wav2 = wmax;
      spec[i].wght = my["g-points"][i].as<Real>();
    }
  } else {
    throw NotFoundError("CorrelatedKTableSpectralGrid",
                        "'g-points' must be defined");
  }
}

SpectralGridPtr SpectralGridFactory::CreateFrom(YAML::Node const& my) {
  if (!my["grid-type"]) {
    throw NotFoundError("SpectralGridFactory", "grid-type");
  }

  std::string grid_type = my["grid-type"].as<std::string>();

  SpectralGridPtr pgrid;

  if (grid_type == "regular") {
    pgrid = std::make_shared<RegularSpacingSpectralGrid>(my);
  } else if (grid_type == "custom") {
    pgrid = std::make_shared<CustomSpacingSpectralGrid>(my);
  } else if (grid_type == "correlated-k-table") {
    pgrid = std::make_shared<CorrelatedKTableSpectralGrid>(my);
  } else {
    throw RuntimeError("SpectralGridFactory", "unknown grid type");
  }

  return pgrid;
}
