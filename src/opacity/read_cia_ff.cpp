// C/C++
#include <exception>
#include <fstream>

// athena
#include <athena/athena.hpp>

// opacity
#include "read_cia_ff.hpp"

std::tuple<AthenaArray<double>, std::vector<double>, std::vector<double>>
read_cia_reform(std::string filename) {
  std::ifstream file{filename};  // open file
  AthenaArray<double> data;      // create storage array
  std::vector<double> temperature_axis;
  std::vector<double> spectral_axis;
  if (file.good()) {
    int nx;              // number of spectral points, horizontal
    int ny;              // number of temperature points, vertical
    double temperature;  // temprary storage for temperature data
    double spectral;     // temprary storage for spectral data
    file >> ny >> nx;
    for (int i = 0; i < ny; ++i) {
      file >> temperature;
      temperature_axis.push_back(temperature);
    }
    for (int i = 0; i < nx; ++i) {
      file >> spectral;
      spectral_axis.push_back(spectral);
    }
    data.NewAthenaArray(ny, nx);
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        file >> data(j, i);
      }
    }
  } else {
    throw std::runtime_error("Unable to open " + filename);
  }
  std::tuple<AthenaArray<double>, std::vector<double>, std::vector<double>>
      data_table = {data, temperature_axis, spectral_axis};
  return data_table;
}

// we are going to read the file twice, first to count # of rows and columns and
// second time
std::tuple<AthenaArray<double>, std::vector<double>, std::vector<double>>
read_freefree(std::string filename) {
  int num_of_row = 0;
  int num_of_column = 1;
  AthenaArray<double> data;  // create storage array
  std::vector<double> temperature_axis;
  std::vector<double> spectral_axis;
  double spectral;
  double temperature;
  std::ifstream cia_file{filename};  // open file
  if (cia_file.good()) {
    std::string line;
    std::getline(cia_file, line);  // get the first line
    std::getline(cia_file, line);  // get the second line
    // calculate the # of columns base on the # of space
    for (int j = 0; j < line.size(); ++j) {
      if (line[j] == ' ') {
        ++num_of_column;
      }
    }
    while (std::getline(cia_file, line)) {
      ++num_of_row;
    }
  } else {
    throw std::runtime_error("Unable to open " + filename);
  }
  cia_file.close();         // close it
  cia_file.open(filename);  // open it again
  if (cia_file.good()) {
    std::string line1;              // storage for the first line
    std::getline(cia_file, line1);  // Skip the first line
    int ny = num_of_row;
    int nx = num_of_column;
    for (int j = 0; j < nx; ++j) {
      cia_file >> temperature;
      temperature_axis.push_back(temperature);  // read off temperature axis
    }
    data.NewAthenaArray(ny, nx);
    for (int j = 0; j < ny; ++j) {
      cia_file >> spectral;  // skip first double and store it into spectal axis
      spectral_axis.push_back(spectral);
      for (int i = 0; i < nx; ++i) {
        cia_file >> data(j, i);
      }
    }
  } else {
    throw std::runtime_error("Unable to open " + filename);
  }
  cia_file.close();
  std::tuple<AthenaArray<double>, std::vector<double>, std::vector<double>>
      data_table = {data, temperature_axis, spectral_axis};
  return data_table;
}
