#ifndef SRC_UTILS_PARAMETER_MAP_HPP_
#define SRC_UTILS_PARAMETER_MAP_HPP_

// C/C++
#include <map>
#include <string>

// external
#include <yaml-cpp/yaml.h>

// athena
#include <athena/athena.hpp>

using ParameterMap = std::unordered_map<std::string, Real>;

inline std::unordered_map<std::string, Real> ToParameterMap(
    const YAML::Node &node) {
  std::unordered_map<std::string, Real> map;
  for (auto it = node.begin(); it != node.end(); ++it) {
    map[it->first.as<std::string>()] = it->second.as<Real>();
  }
  return map;
}

#endif  // SRC_UTILS_PARAMETER_MAP_HPP_
