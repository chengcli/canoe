// C/C++
#include <sstream>
#include <stdexcept>

// canoe
#include <configure.hpp>

// athena
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// climath headers
#include <climath/interpolation.h>

// snap headers
#include <snap/cell_variables.hpp>

// opacity
#include "mwr_absorbers.hpp"
#include "absorption_functions.hpp"

MwrAbsorberH2S::MwrAbsorberH2S(MeshBlock *pmb, ParameterInput *pin,
    int imol, Real xHe):
  Absorber("mw_H2S"), xHe_(xHe)
{
  imols_ = {imol};
  mixrs_ = {1.};

  std::stringstream msg;
  if ((xHe_ < 0.) || (xHe_ > 1.)) {
    msg << "### FATAL ERROR in MwrAbsorberH2S::MwrAbsorberH2S."
        << std::endl << "Value error in molar mixing ratios";
    throw std::runtime_error(msg.str().c_str());
  }
}

MwrAbsorberH2S::MwrAbsorberH2S(MeshBlock *pmb, ParameterInput *pin,
    Real xHe, Real *xH2S, Real *pres, int np):
  Absorber("mw_H2S"), xHe_(xHe)
{
  std::stringstream msg;
  for (int i = 0; i < np; ++i) {
    if ((xH2S[i] < 0.) || (xH2S[i] > 1.)) {
      msg << "### FATAL ERROR in MwrAbsorberNH3::MwrAbsorberNH3."
          << std::endl << "Value error in molar mixing ratios";
      throw std::runtime_error(msg.str().c_str());
    }
    ref_xh2s_.push_back(xH2S[i]);
    ref_pres_.push_back(pres[i]);
  }
}

Real MwrAbsorberH2S::GetAttenuation(Real wave1, Real wave2,
    CellVariables const& var) const
{
  // adapted by cli (Cheng Li), Aug 30
  Real P = var.w[IPR]/1.E5;
  Real T = var.w[IDN];
  Real xdry = 1.;
  for (int i = 1; i <= NVAPOR; ++i) xdry -= var.w[i];
  Real XHe = xHe_*xdry;
  Real XH2, XH2S;

  if (method_ == 1) {
    XH2S = var.w[imols_[0]];
    XH2 = xdry - XHe;
  } else {  // method_ == 2
    XH2S = interp1(var.w[IPR], ref_xh2s_.data(), ref_pres_.data(), ref_pres_.size())*xdry;;
    XH2 = xdry - XHe - XH2S;
  }

  // 1/cm -> 1/m
  return 100.*absorption_coefficient_H2S(wave1, P, T, XH2, XHe, XH2S);
}
