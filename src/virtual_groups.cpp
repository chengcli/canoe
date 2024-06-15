// cantera
#include <cantera/base/yaml.h>

// canoe
#include "virtual_groups.hpp"

void ParameterGroup::SetRealsFrom(YAML::Node const &node) {
  for (auto it = node.begin(); it != node.end(); ++it) {
    params_real_[it->first.as<std::string>()] = it->second.as<Real>();
  }
}

void ParameterGroup::SetIntsFrom(YAML::Node const &node) {
  for (auto it = node.begin(); it != node.end(); ++it) {
    params_int_[it->first.as<std::string>()] = it->second.as<int>();
  }
}
