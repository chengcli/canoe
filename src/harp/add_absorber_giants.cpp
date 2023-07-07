// C/C++
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// external
#include <yaml-cpp/yaml.h>

// application
#include <application/exceptions.hpp>

// opacity
#include <opacity/Giants/microwave/mwr_absorbers.hpp>

// utils
#include <utils/parameter_map.hpp>

// harp
#include "absorber.hpp"
#include "radiation_band.hpp"

namespace gp = GiantPlanets;

std::string extract_name(std::string band_name) {
  // extract name without band type
  std::string delimiter = "-";
  size_t delimiter_pos = aname.find("-");
  return aname.substr(delimiter_pos + delimiter.length());
}

void RadiationBand::addAbsorberGiants(ParameterInput *pin, YAML::Node &node) {
  if (type_ == "radio") {
    addAbsorberGiantsRADIO(pin, node);
  } else if (type_ == "ir") {
    addAbsorberGaintsVISIR(pin, node);
  } else if (type_ == "vis") {
    addAbsorberGaintsVISIR(pin, node);
  } else if (type_ == "uv") {
    addAbsorberGaintsUV(pin, node);
  } else {
    throw NotFoundError("addAbsorberGiants", "Band " + type_);
  }

  absorbers_.back()->SetModel(node["model"].as<std::string>());
}

void RadiationBand::addAbsorberGaintsRADIO(ParameterInput *pin,
                                           YAML::Node &node) {
  auto aname = extract_name(node["name"].as<std::string>());

  if (!node["parameters"]) {
    node["parameters"] = {};
  }

  if (!node["model"]) {
    node["model"] = "unspecified";
  }

  std::vector<std::string> species;

  if (node["dependent-species"])
    species = node["dependent-species"].as<std::vector<std::string>>();

  if (name == "NH3") {
    auto ab = std::make_unique<gp::MwrAbsorberNH3>(
        species, ToParameterMap(node["parameters"]));

    absorbers_.push_back(std::move(ab));
  } else if (name == "H2O") {
    auto ab = std::make_unique<gp::MwrAbsorberH2O>(
        species, ToParameterMap(node["parameters"]));

    absorbers_.push_back(std::move(ab));
  } else if (name == "H2S") {
    auto ab = std::make_unique<gp::MwrAbsorberH2S>(
        species, ToParameterMap(node["parameters"]));

    absorbers_.push_back(std::move(ab));
  } else if (name == "PH3") {
    auto ab = std::make_unique<gp::MwrAbsorberPH3>(
        species, ToParameterMap(node["parameters"]));

    absorbers_.push_back(std::move(ab));
  } else if (name == "CIA") {
    auto ab = std::make_unique<gp::MwrAbsorberCIA>(
        species, ToParameterMap(node["parameters"]));

    absorbers_.push_back(std::move(ab));
  } else if (name == "Electron") {
    auto ab = std::make_unique<gp::MwrAbsorberElectron>(
        species, ToParameterMap(node["parameters"]));

    absorbers_.push_back(std::move(ab));
  } else {
    throw NotFoundError("addAbsorberGiants", "Absorber " + name);
  }
}

void RadiationBand::addAbsorberGiantsVISIR(ParameterInput *pin,
                                           YAML::Node &node) {
  auto name = extract_name(node["name"].as<std::string>());

  if (name == "H2-H2") {
    pabs->addAbsorber(XizH2H2CIA(this, 0, xH2))->loadCoefficient(file);
  } else if (name == "H2-He") {
    pabs->addAbsorber(XizH2HeCIA(this, 0, xH2, xHe))->loadCoefficient(file);
  } else if (name == "freedman_simple") {
    pabs->addAbsorber(FreedmanSimple(this, pin));
  } else if (name == "freedman_mean") {
    pabs->addAbsorber(FreedmanMean(this));
  } else if (strncmp(name.c_str(), "ck-", 3) == 0) {
    auto aname = extract_name(name);
    pabs->addAbsorber(CorrelatedKAbsorber(this, aname));
  } else {
    throw NotFoundError("addAbsorberGiants", "Absorber " + aname);
  }
}

void RadiationBand::addAbsorberGiantsVISIR(ParameterInput *pin,
                                           YAML::Node &node) {
  throw NotImplementedError("addAbsorberGiants", "VISIR");
}
