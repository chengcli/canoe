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
#include <opacity/Giants/freedman.hpp>
#include <opacity/Giants/hydrogen_cia.hpp>
#include <opacity/hitran_absorber.hpp>
#include <opacity/correlatedk_absorber.hpp>

// utils
#include <utils/parameter_map.hpp>

// harp
#include "absorber.hpp"
#include "radiation_band.hpp"

namespace gp = GiantPlanets;

// extract name after hiphen
std::string extract_second(std::string first_second) {
  std::string delimiter = "-";
  size_t delimiter_pos = first_second.find("-");
  return first_second.substr(delimiter_pos + delimiter.length());
}

void RadiationBand::addAbsorberGiants(ParameterInput *pin, YAML::Node &node) {
  if (!node["parameters"]) {
    node["parameters"] = {};
  }

  if (!node["model"]) {
    node["model"] = "unspecified";
  }

  if (type_ == "radio") {
    addAbsorberGiantsRADIO(node);
  } else if (type_ == "ir") {
    addAbsorberGiantsIR(node);
  } else if (type_ == "vis") {
    addAbsorberGiantsVIS(node);
  } else if (type_ == "uv") {
    addAbsorberGiantsUV(node);
  } else {
    throw NotFoundError("addAbsorberGiants", "Band " + type_);
  }

  absorbers_.back()->SetModel(node["model"].as<std::string>());
}

void RadiationBand::addAbsorberGiantsRADIO(YAML::Node &node) {
  auto name = extract_second(node["name"].as<std::string>());
  auto params = ToParameterMap(node["parameters"]);

  std::vector<std::string> species;
  if (node["dependent-species"])
    species = node["dependent-species"].as<std::vector<std::string>>();

  if (name == "NH3") {
    auto ab = std::make_unique<gp::MwrAbsorberNH3>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "H2O") {
    auto ab = std::make_unique<gp::MwrAbsorberH2O>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "H2S") {
    auto ab = std::make_unique<gp::MwrAbsorberH2S>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "PH3") {
    auto ab = std::make_unique<gp::MwrAbsorberPH3>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "CIA") {
    auto ab = std::make_unique<gp::MwrAbsorberCIA>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "Electron") {
    auto ab = std::make_unique<gp::MwrAbsorberElectron>(species, params);

    absorbers_.push_back(std::move(ab));
  } else {
    throw NotFoundError("addAbsorberGiantsRADIO", "Absorber " + name);
  }
}

void RadiationBand::addAbsorberGiantsIR(YAML::Node &node) {
  auto name = extract_second(node["name"].as<std::string>());
  auto params = ToParameterMap(node["parameters"]);

  std::vector<std::string> species;
  if (node["dependent-species"])
    species = node["dependent-species"].as<std::vector<std::string>>();

  if (name == "H2-H2-CIA") {
    auto ab = std::make_unique<XizH2H2CIA>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "H2-He-CIA") {
    auto ab = std::make_unique<XizH2HeCIA>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "CH4") {
    auto ab = std::make_unique<HitranAbsorber>(name, species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "C2H2") {
    auto ab = std::make_unique<HitranAbsorber>(name, species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "C2H4") {
    auto ab = std::make_unique<HitranAbsorber>(name, species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "C2H6") {
    auto ab = std::make_unique<HitranAbsorber>(name, species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "freedman_simple") {
    auto ab = std::make_unique<FreedmanSimple>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "freedman_mean") {
    auto ab = std::make_unique<FreedmanSimple>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (strncmp(name.c_str(), "ck-", 3) == 0) {
    auto aname = extract_second(name);
    auto ab = std::make_unique<CorrelatedKAbsorber>(aname, species, params);

    absorbers_.push_back(std::move(ab));
  } else {
    throw NotFoundError("addAbsorberGiantsIR", "Absorber " + name);
  }
}

void RadiationBand::addAbsorberGiantsVIS(YAML::Node &node) {
  throw NotImplementedError("addAbsorberGiantsVIS");
}

void RadiationBand::addAbsorberGiantsUV(YAML::Node &node) {
  throw NotImplementedError("addAbsorberGiantsUV");
}
