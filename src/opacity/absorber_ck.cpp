// C/C++
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>

// canoe
#include <air_parcel.hpp>
#include <configure.hpp>
#include <constants.hpp>

// climath
#include <climath/interpolation.h>

// harp
#include <harp/spectral_grid.hpp>

// athena
#include <athena/athena.hpp>

// snap
#include <snap/thermodynamics/atm_thermodynamics.hpp>

// opacity
#include "absorber_ck.hpp"

// netcdf
#ifdef NETCDFOUTPUT
extern "C" {
#include <netcdf.h>
}
#endif

void AbsorberCK::LoadCoefficient(std::string fname, int) {
#ifdef NETCDFOUTPUT
    int fileid, varid, err;
    int dimids[10]; // Array to hold dimension IDs (assuming a maximum of 10 dimensions)

    size_t len_p, len_t, len_bin_centers, len_weights;

std::map<std::string, std::string> var_names = {
    {"p", "p"},                        // Pressure (logarithmic, in Pa)
    {"t", "t"},                        // Temperature (K)
    {"bin_centers", "bin_centers"},    // Wavenumber bin centers (cm⁻¹)
    {"bin_edges", "bin_edges"},        // Wavenumber bin edges (cm⁻¹)
    {"weights", "weights"},            // Gaussian quadrature weights
    {"gaussian_points", "gaussian_points"}, // Gaussian quadrature points
    {"kcoeff", "kcoeff"}               // Logarithmic absorption coefficients (log(m²/molecule))
};

    nc_open(fname.c_str(), NC_NOWRITE, &fileid);

    for (const auto& [role, var_name] : var_names) {
        if (role == "kcoeff" || role == "bin_edges") continue; // Handle separately

        nc_inq_varid(fileid, var_name.c_str(), &varid);
        nc_inq_vardimid(fileid, varid, dimids);

        // Assuming each variable is 1D, get the length of the first dimension
        size_t dim_len;
        nc_inq_dimlen(fileid, dimids[0], &dim_len);

        // Assign to len_ based on role
        if (role == "p") {
            len_p = dim_len;
            len_[0] = len_p;
        }
        if (role == "t") {
            len_t = dim_len;
            len_[1] = len_t;
        }
        if (role == "bin_centers") {
            len_bin_centers = dim_len;
            len_[2] = len_bin_centers;
        }
        if (role == "weights") {
            len_weights = dim_len;
            len_[3] = len_weights;
        }
    }

    p_.resize(len_[0]); // Raw pressure
    t_.resize(len_[1]); // Raw temperature
    axis_.resize(len_[0] + len_[1]);             // Pressure + Temperature
    bin_centers_.resize(len_[2]);                 // Bin centers
    weights_.resize(len_[3]);                     // Weights
    bin_edges_.resize(len_[2] + 1);              // Bin edges (len_bin_centers + 1)
    kcoeff_.resize(len_[0] * len_[1] * len_[2] * len_[3]); // Absorption coefficients

    for (const auto& [role, var_name] : var_names) {
        if (role == "kcoeff" || role == "bin_edges") continue; // Handle separately

        // Get the variable ID
        nc_inq_varid(fileid, var_name.c_str(), &varid);

        // Read the data into the corresponding container
        if (role == "p") {
	   nc_get_var_double(fileid, varid, axis_.data());
	   nc_get_var_double(fileid, varid, p_.data());
        }

        if (role == "t") {
            nc_get_var_double(fileid, varid, axis_.data() + len_p);
            nc_get_var_double(fileid, varid, t_.data());
        }
        if (role == "bin_centers") {
            nc_get_var_double(fileid, varid, bin_centers_.data());
        }
        if (role == "weights") {
            nc_get_var_double(fileid, varid, weights_.data());
        }
    }

    nc_inq_varid(fileid, var_names["bin_edges"].c_str(), &varid);
    nc_get_var_double(fileid, varid, bin_edges_.data());

    nc_inq_varid(fileid, var_names["kcoeff"].c_str(), &varid);
    nc_get_var_double(fileid, varid, kcoeff_.data());

    nc_close(fileid);

    // Optional: You can log some debug output to check that everything was loaded correctly
    std::cout << "Loaded coefficient grid: "
              << "Pressure levels = " << len_[0] << ", "
              << "Temperature levels = " << len_[1] << ", "
              << "Wavelength bins = " << len_[2] << ", "
              << "Weights = " << len_[3] << std::endl;
#endif
}

Real AbsorberCK::GetAttenuation(int m, AirParcel const& var) const {
    int iw = m / len_[3]; // Wavelength bin index
    int ig = m % len_[3]; // Weight index
    Real val;
    Real coord[2] = {log(var.w[IPR]), var.w[IDN]};

    // Calculate the starting index for the 2D slice (p, t) at iw, ig
    size_t base_index = iw * len_[3] * len_[0] * len_[1] + ig * len_[0] * len_[1];
    size_t len_2d[2] = {len_[0], len_[1]};

    // Perform 2D interpolation using the flattened kcoeff_
    interpn(&val, coord, &kcoeff_[base_index], axis_.data(), len_2d, 2, 1);

    auto pthermo = Thermodynamics::GetInstance();
    Real dens = var.w[IPR] / (pthermo->GetRd() * var.w[IDN]);

    return exp(val) * dens; // ln(m^2/kg) -> 1/m
}

void AbsorberCK::ModifySpectralGrid(std::vector<SpectralBin>& spec) const {
    // Resize the spec grid to account for both wavelength bins and weights (Gauss points)
    spec.resize(len_[2] * len_[3]); // len_[2] = nbin (wavelength bins), len_[3] = nwght (Gauss points)

    // Loop through both wavelength bins and weights
    for (size_t iw = 0; iw < len_[2]; ++iw) {      // len_[2] is the number of wavelength bins
        for (size_t ig = 0; ig < len_[3]; ++ig) {  // len_[3] is the number of weights (Gauss points)
            size_t m = iw * len_[3] + ig;          // Flattened index for the spectral grid

            // Use bin_edges for wav1 and wav2
            spec[m].wav1 = bin_edges_[iw];       // Lower boundary (wav1)
            spec[m].wav2 = bin_edges_[iw + 1];   // Upper boundary (wav2)
            spec[m].wght = weights_[ig];         // Gauss point weight (from `weights`)
        }
    }
}
