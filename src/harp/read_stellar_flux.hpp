// C/C++
#include <string>
#include <utility>
#include <vector>

//! Read two textfiles and combine the list of data into a pair of vector of
//! double.
//! Example files: sw_band_flux_HD189_11.txt and wavelengths_GCM_11.txt
std::pair<std::vector<double>, std::vector<double>> read_stellar_flux(
    std::string file1, std::string file2);
