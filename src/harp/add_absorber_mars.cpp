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
#include <opacity/Giants/freedman.hpp>
#include <opacity/Giants/hydrogen_cia.hpp>
#include <opacity/Giants/microwave/mwr_absorbers.hpp>
#include <opacity/correlatedk_absorber.hpp>
#include <opacity/hitran_absorber.hpp>

// utils
#include <utils/parameter_map.hpp>

// harp
#include "absorber.hpp"
#include "radiation_band.hpp"

void RadiationBand::addAbsorberMars(ParameterInput *pin, YAML::Node &node) {
  if (!node["parameters"]) {
    node["parameters"] = {};
  }

  if (!node["model"]) {
    node["model"] = "unspecified";
  }

  if (category_ == "infrared") {
    addAbsorberMarsInfrared(node);
  } else {
    throw NotFoundError("addAbsorberMars", "Category: " + category_);
  }

  absorbers_.back()->SetModel(node["model"].as<std::string>());
}

void RadiationBand::addAbsorberMarsInfrared(YAML::Node &node) {
  auto name = node["name"].as<std::string>();
  auto params = ToParameterMap(node["parameters"]);

  std::vector<std::string> species;
  if (node["dependent-species"])
    species = node["dependent-species"].as<std::vector<std::string>>();

  if (name == "CO2-CO2-CIA") {
    // auto ab = std::make_unique<XizH2H2CIA>(species, params);

    // absorbers_.push_back(std::move(ab));
  } else if (name == "CO2") {
    auto ab = std::make_unique<HitranAbsorber>(name, species, params);
    absorbers_.push_back(std::move(ab));
  } else if (name == "H2O") {
    auto ab = std::make_unique<HitranAbsorber>(name, species, params);
    absorbers_.push_back(std::move(ab));
  } else if (name == "O3") {
    auto ab = std::make_unique<HitranAbsorber>(name, species, params);
    absorbers_.push_back(std::move(ab));
  } else if (name == "O2") {
    auto ab = std::make_unique<HitranAbsorber>(name, species, params);
    absorbers_.push_back(std::move(ab));
  } else if (name == "CO") {
    auto ab = std::make_unique<HitranAbsorber>(name, species, params);
    absorbers_.push_back(std::move(ab));
  } else {
    throw NotFoundError("addAbsorberMarsInfrared", "Absorber " + name);
  }
}
