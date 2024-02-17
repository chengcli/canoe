// C/C++
#include <fstream>
#include <utility>
#include <exception>

// opacity
#include "read_rayleigh.hpp"

std::vector<double> read_rayleigh(std::string file) {
  std::ifstream data{file};                  // open file
  std::vector<double> cross_section_output;  // store the output
  std::string line;  // a temporary storage to get over the first line
  // if stream OK = file readable
  if (data.good()) {
    double x;
    std::getline(data, line);  // Skip the first line
    // as long as next value readable
    while (data >> x) {
      // put it into output
      cross_section_output.push_back(x);
    }
  } else {
    throw std::runtime_error("Unable to open " + file);
  }
  return cross_section_output;
}
