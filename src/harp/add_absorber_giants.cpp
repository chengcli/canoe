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

void RadiationBand::addAbsorberGiants(ParameterInput *pin, YAML::Node &node) {
  std::string aname = node["name"].as<std::string>();

  // extract name without band type
  std::string delimiter = "-";
  size_t delimiter_pos = aname.find("-");
  aname = aname.substr(delimiter_pos + delimiter.length());

  if (!node["parameters"]) {
    node["parameters"] = {};
  }

  if (!node["model"]) {
    node["model"] = "unspecified";
  }

  std::vector<std::string> species;
  if (node["dependent-species"])
    species = node["dependent-species"].as<std::vector<std::string>>();

  if (type_ == "radio") {
    if (aname == "NH3") {
      auto ab = std::make_unique<gp::MwrAbsorberNH3>(
          pmy_block_, species, ToParameterMap(node["parameters"]));

      absorbers_.push_back(std::move(ab));
    } else if (aname == "H2O") {
      auto ab = std::make_unique<gp::MwrAbsorberH2O>(
          pmy_block_, species, ToParameterMap(node["parameters"]));

      absorbers_.push_back(std::move(ab));
    } else if (aname == "H2S") {
      auto ab = std::make_unique<gp::MwrAbsorberH2S>(
          pmy_block_, species, ToParameterMap(node["parameters"]));

      absorbers_.push_back(std::move(ab));
    } else if (aname == "PH3") {
      auto ab = std::make_unique<gp::MwrAbsorberPH3>(
          pmy_block_, species, ToParameterMap(node["parameters"]));

      absorbers_.push_back(std::move(ab));
    } else if (aname == "CIA") {
      auto ab = std::make_unique<gp::MwrAbsorberCIA>(
          pmy_block_, species, ToParameterMap(node["parameters"]));

      absorbers_.push_back(std::move(ab));
    } else if (aname == "Electron") {
      auto ab = std::make_unique<gp::MwrAbsorberElectron>(
          pmy_block_, species, ToParameterMap(node["parameters"]));

      ab->SetModel(node["model"].as<std::string>());
      absorbers_.push_back(std::move(ab));
    } else {
      throw NotFoundError("addAbsorberGiants", "Absorber " + aname);
    }
  } else if (type_ == "ir") {
    auto ab = std::make_unique<Absorber>(aname);
    absorbers_.push_back(std::move(ab));
  } else if (type_ == "vis") {
    auto ab = std::make_unique<Absorber>(aname);
    absorbers_.push_back(std::move(ab));
  } else if (type_ == "uv") {
    auto ab = std::make_unique<Absorber>(aname);
    absorbers_.push_back(std::move(ab));
  } else {
    throw NotFoundError("addAbsorberGiants", "Band " + type_);
  }
}
