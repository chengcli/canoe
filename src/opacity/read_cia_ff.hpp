// C/C++
#include <string>
#include <tuple>
#include <utility>
#include <vector>

template <typename T>
class AthenaArray;

//! The first array is the 2D data sheet of cross sections, the second vector is
//! the temperature axis, and the third is spectral axis

//! read cia reform format file on temperature vs. spectral wavelength 2D data
//! set.
std::tuple<AthenaArray<double>, std::vector<double>, std::vector<double>>
read_cia_reform(std::string filename);

//! read free-free absorption format file on temperature vs. spectral wavelength
//! 2D data set.
std::tuple<AthenaArray<double>, std::vector<double>, std::vector<double>>
read_freefree(std::string filename);
