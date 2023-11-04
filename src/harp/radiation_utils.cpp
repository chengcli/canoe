// C++
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

// application
#include <application/exceptions.hpp>

// canoe
#include <configure.hpp>

// climath
#include <climath/core.h>

// harp
#include "radiation_utils.hpp"

/*void Radiation::TotalFlux(AthenaArray<Real>& flux) const
{
  // count how many bands
  int b = 0;
  int max_ntau = 0;
  Radiation const *p = this;
  while (p != NULL) {
    max_ntau = std::max(max_ntau, p->ntau);
    p = p->next;
    b++;
  }

  // check dimension consistancy, reallocate memory if necessary
  if (flux.GetDim1() != b && flux.GetDim2() < max_ntau) {
    flux.DeleteAthenaArray();
    flux.NewAthenaArray(max_ntau, b);
  }

  b = 0;
  p = this;
  while (p != NULL) {
    for (int j = 0; j < p->ntau; ++j)
      flux(j,b) = 0.;
    for (int i = 0; i < p->nwave; ++i)
      for (int j = 0; j < p->ntau; ++j)
        // flup - rfldn - rfldir
        flux(j,b) += p->weight_[i]*(p->rad(i,j,2) - p->rad(i,j,1) -
p->rad(i,j,0)); p = p->next; b++;
  }
}

void Radiation::SumTotalFlux(AthenaArray<Real>& tflux, AthenaArray<Real>&
disort_flux, int gk, int gj) const
{
  // count how many bands
  int b = 0;
  int max_ntau = 0;
  Radiation const *p = this;
  while (p != NULL) {
    max_ntau = std::max(max_ntau, p->ntau);
    p = p->next;
    b++;
  }
  AthenaArray<Real> flux;
  flux.NewAthenaArray(max_ntau, b);
  b = 0;
  p = this;
  while (p != NULL) {
    for (int j = 0; j < p->ntau; ++j)
      flux(j,b) = 0.;
    for (int j = 0; j < p->ntau+2*(NGHOST); ++j){
      disort_flux(0,gk,gj,j) = 0.;
      disort_flux(1,gk,gj,j) = 0.;
      disort_flux(2,gk,gj,j) = 0.;
    }
    for (int i = 0; i < p->nwave; ++i)
      for (int j = 0; j < p->ntau; ++j){
        // flup - rfldn - rfldir
        disort_flux(0,gk,gj,j+NGHOST) += p->weight_[i]*p->rad(i,j,2);
        disort_flux(1,gk,gj,j+NGHOST) += p->weight_[i]*p->rad(i,j,1);
        disort_flux(2,gk,gj,j+NGHOST) += p->weight_[i]*p->rad(i,j,0);
        flux(j,b) += p->weight_[i]*(p->rad(i,j,2) - p->rad(i,j,1) -
p->rad(i,j,0));
      }
    p = p->next;
    b++;
  }

  // sum; thl
  for (int j =0; j < max_ntau; ++j){
      tflux(j) = 0.;
  }
  for (int j =0; j <max_ntau; ++j)
    for (int bb =0; bb <b; ++bb){
      tflux(j) += flux(max_ntau-j-1,bb);
    }
}

void Radiation::WriteTopFlux(std::string fname) const
{
  std::cout << "Top flux written into file: " << fname << std::endl;
  std::ofstream out(fname.c_str(), std::ios::out);
  out << std::left << std::setw(20) << std::setprecision(8) << "#
Wavenumber(cm-1)"
      << std::left << std::setw(20) << "Top rfldir[]"
      << std::left << std::setw(20) << "Top rfldn[]"
      << std::left << std::setw(20) << "Top flup[]"
      << std::endl;

  Radiation const *p = this;
  while (p != NULL) {
    for (int i = 0; i < p->nwave; ++i)
      out << std::left << std::setw(20) << p->wave[i]
          << std::left << std::setw(20) << p->rad(i,0,0)
          << std::left << std::setw(20) << p->rad(i,0,1)
          << std::left << std::setw(20) << p->rad(i,0,2)
          << std::endl;
    p = p->next;
  }
}

void Radiation::WriteOpticalDepth(std::string fname) const
{
  std::cout << "Optical depth written into file: " << fname << std::endl;
  std::ofstream out(fname.c_str(), std::ios::out);

  Radiation const *p = this;
  while (p != NULL) {
    out << std::left << "# Band: " << p->myname << std::endl;
    out << std::setw(40) << "# Number of wavenumbers : " << p->nwave <<
std::endl; out << std::setw(40) << "# Number of levels: " << p->nlevel <<
std::endl; out << std::right; for (int i = 0; i < p->nlevel; ++i) out <<
std::setw(16) << std::setprecision(8) << p->level[i]/1.E3; out << std::endl;

    for (int i = 0; i < p->nwave; ++i) {
      out << std::setw(16) << std::setprecision(8) << p->wave[i];
      for (int j = 0; j < p->nlevel - 1; ++j)
        out << std::setw(16) << std::setprecision(8) << p->tau(ITAU,i,j);
      out << std::endl;
    }
    p = p->next;
  }
}

void Radiation::WriteTopRadiance(std::string fname) const
{
  std::cout << "Top radiance written into file: " << fname << std::endl;
  std::ofstream out(fname.c_str(), std::ios::out);
  out << std::left << std::setw(20) << std::setprecision(8)
      << "# Wavenumber/Frequency" << std::endl;

  Radiation const *p = this;
  while (p != NULL) {
    for (int i = 0; i < p->nwave; ++i) {
      out << std::left << std::setw(20) << p->wave[i];
      for (int j = 0; j < p->nphi; ++j)
        for (int k = 0; k < p->numu; ++k)
          out << std::left << std::setw(20) << p->uu(i,j,0,k);
      out << std::endl;
    }
    p = p->next;
  }
}

void WriteHeatingRate(std::string fname, AthenaArray<Real> const& flux,
      AthenaArray<Real> const& hrate, Real const* level)
{
  int nband = flux.GetDim1();
  int nlevel = flux.GetDim2();

  // print mean intensity, flux divergence and heating rate
  std::cout << "Heating rate written into file: " << fname << std::endl;
  std::ofstream out(fname.c_str(), std::ios::out);
  out  << std::left << std::setw(6) << "Level"
       << std::left << std::setw(12) << "Height [km]";

  for (int b = 0; b < nband; ++b) {
    char bname[80];
    sprintf(bname, "B%d flux [w/m^2]", b + 1);
    out << std::left << std::setw(16) << bname
        << std::left << std::setw(20) << "Heating rate [w/kg]";
  }

  if (nband > 1) {
    out << std::left << std::setw(16) << "Flux [w/m^2]"
        << std::left << std::setw(20) << "Heating rate [w/kg]"
        << std::endl;
  } else {
    out << std::endl;
  }

  for (int i = 0; i < nlevel; ++i) {
    out  << std::left << std::setw(6) << i + 1
         << std::left << std::setw(12) << level[i]/1.E3;

    Real total_flux = 0., total_hrate = 0.;
    for (int b = 0; b < nband; ++b) {
      out << std::left << std::setw(16) << std::setprecision(8) << flux(i,b)
          << std::left << std::setw(20) << std::setprecision(8) << hrate(i,b);
      total_flux += flux(i,b);
      total_hrate += hrate(i,b);
    }

    if (nband > 1) {
      out << std::left << std::setw(16) << std::setprecision(8) << total_flux
          << std::left << std::setw(20) << std::setprecision(8) << total_hrate
          << std::endl;
    } else {
      out << std::endl;
    }
  }
}*/

void packSpectralProperties(Real *buf, Real const *tau, Real const *ssa,
                            Real const *pmom, int nlayer, int npmom) {}

void unpackSpectralProperties(Real *tau, Real *ssa, Real *pmom, Real const *buf,
                              int slyr, int npmom, int nblocks, int npmom_max) {
}
