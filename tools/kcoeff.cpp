// C/C++
#include <unistd.h>

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// canoe
#include <configure.hpp>

#ifdef NETCDFOUTPUT
extern "C" {
#include "netcdf.h"
}
#endif

template <typename T, int N>
class Array {
 protected:
  int m_size[N];
  std::vector<T> m_data;

 public:
  Array() {}

  void resize(int n2, int n1) {
    m_size[0] = n1;
    m_size[1] = n2;

    m_data.resize(n1 * n2);
  }

  void resize(int n3, int n2, int n1) {
    m_size[0] = n1;
    m_size[1] = n2;
    m_size[2] = n3;

    m_data.resize(n1 * n2 * n3);
  }

  T* data() { return m_data.data(); }

  T const* data() const { return m_data.data(); }

  T& operator()(int i2, int i1) { return m_data[i2 * m_size[0] + i1]; }

  T const& operator()(int i2, int i1) const {
    return m_data[i2 * m_size[0] + i1];
  }

  T& operator()(int i3, int i2, int i1) {
    return m_data[i3 * m_size[1] * m_size[0] + i2 * m_size[0] + i1];
  }

  T const& operator()(int i3, int i2, int i1) const {
    return m_data[i3 * m_size[1] * m_size[0] + i2 * m_size[0] + i1];
  }
};

bool fexists(std::string fname) {
  std::ifstream ifile(fname.c_str());
  return ifile.is_open();
}

std::string decomment(std::string fname) {
  if (!fexists(fname)) {
    std::cout << fname << " does not exist." << std::endl;
    assert(0);
  }
  // assert(fexists(fname));
  std::ifstream file(fname.c_str(), std::ios::in);
  std::string ss;
  char c;
  while (file) {
    file.get(c);
    if (c == '#') {
      while (c != '\n' && file) file.get(c);
      continue;
    }
    ss += c;
  }
  return ss;
}

/** @file
 * @brief Read absorption data file and interpolate to a given temperature and
 * wavenumber grid.
 *
 * Molecular absorption coefficients are calculated by RFM
 */

int main(int argc, char* argv[]) {
  clock_t cpu_tic = clock();
  std::cout << "kcoeff program running ..." << std::endl;

  std::string line;
  std::vector<std::string> mols, cfiles, mfiles;
  std::vector<std::vector<double> > cwave, ctemp;
  std::vector<double> wave_axis, temp_axis, temp, pres;
  Array<double, 3> kcoeff_mol;
  int nmols, nwaves, ntemps, nlevels;
  double wmin, wmax, tmin, tmax, junk;

  int iarg;
  std::string inpfile, outfile;
  while ((iarg = getopt(argc, argv, "i:o:")) != -1) {
    switch (iarg) {
      case 'i':
        inpfile = std::string(optarg);
        break;
      case 'o':
        outfile = std::string(optarg);
        break;
      default:
        std::cerr << "unhandled option." << std::endl;
        exit(1);
    }
  }

  std::stringstream inp(decomment(inpfile));

  // molecule absorber
  inp >> nmols;
  mols.resize(nmols);
  mfiles.resize(nmols);
  for (int i = 0; i < nmols; i++) inp >> mols[i];
  for (int i = 0; i < nmols; i++) inp >> mfiles[i];

  // wavenumber
  inp >> wmin >> wmax >> nwaves;
  wave_axis.resize(nwaves);
  for (int i = 0; i < nwaves; i++)
    wave_axis[i] = wmin + i * (wmax - wmin) / (nwaves - 1);

  // relative temperature
  inp >> tmin >> tmax >> ntemps;
  temp_axis.resize(ntemps);
  for (int i = 0; i < ntemps; i++)
    temp_axis[i] = tmin + i * (tmax - tmin) / (ntemps - 1);

  // number of levels
  inp >> nlevels;
  temp.resize(nlevels);
  pres.resize(nlevels);

  // temperature profile
  for (int i = 0; i < nlevels; i++) inp >> temp[i];

  // pressure profile
  for (int i = 0; i < nlevels; i++) {
    inp >> pres[i];
    pres[i] *= 100.;  // mb -> pa
  }

  if (outfile.empty()) {
    std::cerr << "output netcdf file not specifeid." << std::endl;
    exit(1);
  }

#ifdef NETCDFOUTPUT

  // write netcdf output
  int ncid, dimid, varid;
  nc_create(outfile.c_str(), NC_NETCDF4, &ncid);

  int wave_dimid, pres_dimid, temp_dimid;
  int wave_varid, pres_varid, temp_varid;

  nc_def_dim(ncid, "Wavenumber", wave_axis.size(), &wave_dimid);
  nc_def_var(ncid, "Wavenumber", NC_DOUBLE, 1, &wave_dimid, &wave_varid);
  nc_put_att_text(ncid, wave_varid, "units", 4, "1/cm");
  nc_put_att_text(ncid, wave_varid, "long_name", 15, "wavenumber grid");

  nc_def_dim(ncid, "Pressure", pres.size(), &pres_dimid);
  nc_def_var(ncid, "Pressure", NC_DOUBLE, 1, &pres_dimid, &pres_varid);
  nc_put_att_text(ncid, pres_varid, "units", 2, "pa");
  nc_put_att_text(ncid, pres_varid, "long_name", 18, "reference pressure");

  nc_def_dim(ncid, "TempGrid", temp_axis.size(), &temp_dimid);
  nc_def_var(ncid, "TempGrid", NC_DOUBLE, 1, &temp_dimid, &temp_varid);
  nc_put_att_text(ncid, temp_varid, "units", 1, "K");
  nc_put_att_text(ncid, temp_varid, "long_name", 24,
                  "temperature anomaly grid");

  int tp_varid, *mol_varid;
  mol_varid = new int[nmols];
  nc_def_var(ncid, "Temperature", NC_DOUBLE, 1, &pres_dimid, &tp_varid);
  nc_put_att_text(ncid, tp_varid, "long_name", 21, "reference temperature");
  nc_put_att_text(ncid, tp_varid, "units", 1, "K");

  int dims[3] = {wave_dimid, pres_dimid, temp_dimid};
  for (int i = 0; i < nmols; i++) {
    nc_def_var(ncid, mols[i].c_str(), NC_DOUBLE, 3, dims, mol_varid + i);
    nc_put_att_text(ncid, mol_varid[i], "units", 12, "ln(m^2/kmol)");
    nc_put_att_text(ncid, mol_varid[i], "long_name", 24,
                    "absorption cross section");
  }

  nc_enddef(ncid);

  nc_put_var_double(ncid, wave_varid, wave_axis.data());
  nc_put_var_double(ncid, pres_varid, pres.data());
  nc_put_var_double(ncid, tp_varid, temp.data());
  nc_put_var_double(ncid, temp_varid, temp_axis.data());

  std::ifstream infile;

  // molecule absorption
  kcoeff_mol.resize(nwaves, nlevels, ntemps);
  for (int m = 0; m < nmols; m++) {
    std::cout << "making absorption file for " << mols[m] << " ..."
              << std::endl;
    infile.open(mfiles[m].c_str());
    if (infile.is_open()) {
      for (int i = 0; i < 10; i++) getline(infile, line);
      for (int k = 0; k < nwaves; k++) {
        infile >> junk;
        for (int j = 0; j < nlevels; j++)
          for (int i = 0; i < ntemps; i++)
            infile >> kcoeff_mol(k, nlevels - j - 1, i);
      }
      infile.close();
    } else {
      std::cerr << "cannot open file: " << mfiles[m] << std::endl;
      exit(1);
    }
    nc_put_var_double(ncid, mol_varid[m], kcoeff_mol.data());
  }

  nc_close(ncid);

  delete[] mol_varid;

#else   // NO_NETCDFOUTPUT
  std::cout << "NO NETCDF, NO OUPUT" << std::endl;
#endif  // NETCDFOUTPUT

  clock_t cpu_toc = clock();
  std::cout << "kcoeff program ends successfully in "
            << double(cpu_toc - cpu_tic) / CLOCKS_PER_SEC << " s." << std::endl
            << "output file written into " << outfile << std::endl;
}
