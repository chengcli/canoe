// C/C++
#include "read_stellar_flux.hpp"

#include <fstream>

std::pair<std::vector<double>, std::vector<double>> read_stellar_flux(
    std::string file1, std::string file2) {
  std::ifstream file_flux{file1};        // open file
  std::ifstream file_wavelength{file2};  // open file
  std::vector<double> flux_output;
  std::vector<double> wavelength_output;
  // if stream OK = file readable
  if (file_flux.good()) {
    double x;
    // as long as next value readable
    while (file_flux >> x) {
      // put it into output1
      flux_output.push_back(x);
    }
  } else {
    throw std::runtime_error("Unable to open " + file1);
  }

  if (file_wavelength.good()) {
    double y;
    file_wavelength.ignore(2);
    // as long as next value readable
    while (file_wavelength >> y) {
      // put it into output2
      wavelength_output.push_back(y);
    }
  } else {
    throw std::runtime_error("Unable to open " + file2);
  }
  std::pair<std::vector<double>, std::vector<double>> combined_output = {
      flux_output, wavelength_output};
  return combined_output;
}
