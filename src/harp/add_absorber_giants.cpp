// C/C++
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

void RadiationBand::addAbsorberGiants(ParameterInput *pin, std::string bname,
                                      YAML::Node &node) {
  std::string aname = node["name"].as<std::string>();

  if (bname == "radio") {
    std::vector<std::string> species;

    if (aname == "NH3" || aname == "radio-NH3") {
      if (node["dependent-species"])
        species = node["dependent-species"].as<std::vector<std::string>>();

      auto ab = std::make_unique<gp::MwrAbsorberNH3>(
          pmy_block_, species, ToParameterMap(node["parameters"]));
      ab->SetModel(node["model"].as<std::string>());

      absorbers_.push_back(std::move(ab));
    } else if (aname == "H2O" || aname == "radio-H2O") {
      if (node["dependent-species"])
        species = node["dependent-species"].as<std::vector<std::string>>();

      auto ab = std::make_unique<gp::MwrAbsorberH2O>(
          pmy_block_, species, ToParameterMap(node["parameters"]));
      ab->SetModel(node["model"].as<std::string>());

      absorbers_.push_back(std::move(ab));
    } else if (aname == "H2S" || aname == "radio-H2S") {
      if (node["dependent-species"])
        species = node["dependent-species"].as<std::vector<std::string>>();

      auto ab = std::make_unique<gp::MwrAbsorberH2S>(
          pmy_block_, species, ToParameterMap(node["parameters"]));

      absorbers_.push_back(std::move(ab));
    } else if (aname == "PH3" || aname == "radio-PH3") {
      if (node["dependent-species"])
        species = node["dependent-species"].as<std::vector<std::string>>();

      auto ab = std::make_unique<gp::MwrAbsorberPH3>(
          pmy_block_, species, ToParameterMap(node["parameters"]));

      absorbers_.push_back(std::move(ab));
    } else if (aname == "CIA" || aname == "radio-CIA") {
      if (node["dependent-species"])
        species = node["dependent-species"].as<std::vector<std::string>>();

      auto ab = std::make_unique<gp::MwrAbsorberCIA>(
          pmy_block_, species, ToParameterMap(node["parameters"]));

      absorbers_.push_back(std::move(ab));
    } else if (aname == "Electron" || aname == "radio-Electron") {
      if (node["dependent-species"])
        species = node["dependent-species"].as<std::vector<std::string>>();

      auto ab = std::make_unique<gp::MwrAbsorberElectron>(
          pmy_block_, species, ToParameterMap(node["parameters"]));

      absorbers_.push_back(std::move(ab));
    } else {
      throw NotFoundError("addAbsorberGiants", "Absorber " + aname);
    }
  } else if (bname == "ir") {
    auto ab = std::make_unique<Absorber>(aname);
    absorbers_.push_back(std::move(ab));
  } else if (bname == "vis") {
    auto ab = std::make_unique<Absorber>(aname);
    absorbers_.push_back(std::move(ab));
  } else if (bname == "uv") {
    auto ab = std::make_unique<Absorber>(aname);
    absorbers_.push_back(std::move(ab));
  } else {
    throw NotFoundError("addAbsorberGiants", "Band " + bname);
  }
}
