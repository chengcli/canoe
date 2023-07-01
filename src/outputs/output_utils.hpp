#ifndef SRC_OUTPUTS_OUTPUT_UTILS_HPP_
#define SRC_OUTPUTS_OUTPUT_UTILS_HPP_

// C/C++
#include <string>
#include <vector>

// athena
#include <athena/athena.hpp>

int get_num_variables(std::string grid, AthenaArray<Real> const& data);

std::string get_grid_type(std::string name);

std::string get_units(std::string name);

std::string get_long_name(std::string name);

using DiagnosticTable = std::vector<std::vector<std::string>>;

#endif  // SRC_OUTPUTS_OUTPUT_UTILS_HPP_
