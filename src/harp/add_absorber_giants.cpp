// C/C++
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// external
#include <yaml-cpp/yaml.h>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// opacity
#include <opacity/Giants/freedman.hpp>
#include <opacity/Giants/hydrogen_cia.hpp>
#include <opacity/Giants/microwave/mwr_absorbers.hpp>
#include <opacity/correlatedk_absorber.hpp>
#include <opacity/hitran_absorber.hpp>

// utils
#include <utils/extract_substring.hpp>
#include <utils/fileio.hpp>
#include <utils/parameter_map.hpp>

// harp
#include "absorber.hpp"
#include "radiation_band.hpp"

namespace gp = GiantPlanets;

void RadiationBand::addAbsorberGiants(ParameterInput *pin, YAML::Node &node) {
  if (!node["parameters"]) {
    node["parameters"] = {};
  }

  if (!node["model"]) {
    node["model"] = "unspecified";
  }

  if (category_ == "radio") {
    addAbsorberGiantsRadio(node);
  } else if (category_ == "infrared") {
    addAbsorberGiantsInfrared(node);
  } else if (category_ == "visible") {
    addAbsorberGiantsVisible(node);
  } else if (category_ == "ultraviolet") {
    addAbsorberGiantsUltraviolet(node);
  } else {
    throw NotFoundError("addAbsorberGiants", "Category: " + category_);
  }

  absorbers_.back()->SetModel(node["model"].as<std::string>());
}

void RadiationBand::addAbsorberGiantsRadio(YAML::Node &node) {
  auto name = node["name"].as<std::string>();
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
    throw NotFoundError("addAbsorberGiantsRadio", "Absorber " + name);
  }
}

void RadiationBand::addAbsorberGiantsInfrared(YAML::Node &node) {
  auto name = node["name"].as<std::string>();
  auto params = ToParameterMap(node["parameters"]);

  std::vector<std::string> species;
  std::string data_file;

  if (node["dependent-species"])
    species = node["dependent-species"].as<std::vector<std::string>>();

  if (node["data"]) {
    data_file = node["data"].as<std::string>();
  } else {
    data_file = "kcoeff." + myfile_ + "-" + name_ + ".nc";
  }

  auto app = Application::GetInstance();
  auto full_path = app->FindInputFile(data_file);

  if (name == "H2-H2-CIA") {
    auto ab = std::make_unique<XizH2H2CIA>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "H2-He-CIA") {
    auto ab = std::make_unique<XizH2HeCIA>(species, params);

    absorbers_.push_back(std::move(ab));
  } else if (name == "CH4") {
    auto ab = std::make_unique<HitranAbsorber>(name, species, params);
    if (FileExists(full_path)) ab->LoadCoefficient(data_file);
    absorbers_.push_back(std::move(ab));

  } else if (name == "C2H2") {
    auto ab = std::make_unique<HitranAbsorber>(name, species, params);
    if (FileExists(full_path)) ab->LoadCoefficient(data_file);
    absorbers_.push_back(std::move(ab));

  } else if (name == "C2H4") {
    auto ab = std::make_unique<HitranAbsorber>(name, species, params);
    if (FileExists(full_path)) ab->LoadCoefficient(data_file);
    absorbers_.push_back(std::move(ab));

  } else if (name == "C2H6") {
    auto ab = std::make_unique<HitranAbsorber>(name, species, params);
    if (FileExists(full_path)) ab->LoadCoefficient(data_file);
    absorbers_.push_back(std::move(ab));

  } else if (name == "freedman_simple") {
    auto ab = std::make_unique<FreedmanSimple>(species, params);
    absorbers_.push_back(std::move(ab));

  } else if (name == "freedman_mean") {
    auto ab = std::make_unique<FreedmanSimple>(species, params);
    absorbers_.push_back(std::move(ab));

  } else if (strncmp(name.c_str(), "ck-", 3) == 0) {
    auto aname = extract_second(name, "-");
    auto ab = std::make_unique<CorrelatedKAbsorber>(aname, species, params);

    absorbers_.push_back(std::move(ab));
  } else {
    throw NotFoundError("addAbsorberGiantsIR", "Absorber " + name);
  }
}

void RadiationBand::addAbsorberGiantsVisible(YAML::Node &node) {
  throw NotImplementedError("addAbsorberGiantsVisible");
}

void RadiationBand::addAbsorberGiantsUltraviolet(YAML::Node &node) {
  throw NotImplementedError("addAbsorberGiantsUV");
}
