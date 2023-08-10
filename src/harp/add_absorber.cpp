// external
#include <yaml-cpp/yaml.h>

// athena++
#include <athena/parameter_input.hpp>

// harp
#include "radiation_band.hpp"

// all bands
void RadiationBand::AddAbsorber(ParameterInput *pin, YAML::Node &node) {
  if (!node["parameters"]) {
    node["parameters"] = {};
  }

  if (!node["model"]) {
    node["model"] = "unspecified";
  }

  if (category_ == "radio") {
    addAbsorberRadio(node);
  } else if (category_ == "infrared") {
    addAbsorberInfrared(node);
  } else if (category_ == "visible") {
    addAbsorberVisible(node);
  } else if (category_ == "ultraviolet") {
    addAbsorberUltraviolet(node);
  } else {
    throw NotFoundError("addAbsorberGiants", "Category: " + category_);
  }

  absorbers_.back()->SetModel(node["model"].as<std::string>());
}

// radio bands
void __attribute__((weak)) RadiationBand::addAbsorberRadio(YAML::Node &node) {}

// infrared bands
void __attribute__((weak))
RadiationBand::addAbsorberInfrared(YAML::Node &node) {}

// visible bands
void __attribute__((weak)) RadiationBand::addAbsorberVisible(YAML::Node &node) {
}

// uv bands
void __attribute__((weak))
RadiationBand::addAbsorberUltraviolet(YAML::Node &node) {}
