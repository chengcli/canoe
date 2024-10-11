// C/C++
#include <fstream>
#include <iostream>
#include <string>

// athena
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// external
#include <yaml-cpp/yaml.h>

#include <application/exceptions.hpp>

// microphysics
#include "microphysical_schemes.hpp"
#include "microphysics.hpp"

AllMicrophysicalSchemes MicrophysicalSchemesFactory::Create(
    MeshBlock *pmb, ParameterInput *pin) {
  std::string input_key = Microphysics::input_key;

  AllMicrophysicalSchemes systems;

  if (pin->DoesParameterExist("chemistry", input_key)) {
    std::string filename = pin->GetString("chemistry", input_key);
    std::ifstream stream(filename);
    if (stream.good() == false) {
      throw NotFoundError("microphysics",
                          "Cannot open microphysics config file: " + filename);
    }

    YAML::Node node = YAML::Load(stream);
    if (!node["microphysics"]) {
      throw NotFoundError("Microphysics", "microphysics");
    }

    for (auto sys : node["microphysics"]) {
      std::string name = sys.as<std::string>();
      std::string scheme = node[name]["scheme"].as<std::string>();
      if (scheme == "Kessler94") {
        // FIXME(cli)
        // auto p = std::make_shared<Kessler94>(name, node[name]);
        // systems.push_back(p);
      } else {
        throw NotFoundError("Microphysics", scheme);
      }
    }
  }

  return systems;
}
