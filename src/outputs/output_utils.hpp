#ifndef OUTPUT_UTILS_HPP
#define OUTPUT_UTILS_HPP

// C/C++ headers
#include <string>

// Athena++ headers
#include <athena/athena.hpp>

int get_num_variables(std::string grid, AthenaArray<Real> const& data);

std::string get_grid_type(std::string name);

std::string get_units(std::string name);

std::string get_long_name(std::string name);

#endif
