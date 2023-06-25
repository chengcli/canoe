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

// snap
#include <snap/cell_variables.hpp>

// opacity
#include "mwr_absorbers.hpp"
#include "absorption_functions.hpp"

MwrAbsorberNH3::MwrAbsorberNH3(MeshBlock *pmb, ParameterInput *pin,
    int imol, Real xHe, Real xH2O):
  Absorber("mw_NH3"), method_(1), xHe_(xHe), xH2O_(xH2O)
{
  imols_ = {imol};
  mixrs_ = {1.};

  std::stringstream msg;
  if ((xHe < 0.) || (xH2O < 0.) || (xHe + xH2O > 1.)) {
    msg << "### FATAL ERROR in MwrAbsorberNH3::MwrAbsorberNH3."
        << std::endl << "Value error in molar mixing ratios";
    ATHENA_ERROR(msg);
  }
}

MwrAbsorberNH3::MwrAbsorberNH3(MeshBlock *pmb, ParameterInput *pin,
    int imol, Real xHe, Real *xH2O, Real *pres, int np):
  Absorber("mw_NH3"), method_(2), xHe_(xHe)
{
  imols_ = {imol};
  mixrs_ = {1.};

  std::stringstream msg;
  for (int i = 0; i < np; ++i) {
    if ((xH2O[i] < 0.) || (xH2O[i] > 1.)) {
      msg << "### FATAL ERROR in MwrAbsorberNH3::MwrAbsorberNH3."
          << std::endl << "Value error in molar mixing ratios";
      ATHENA_ERROR(msg);
    }
    ref_xh2o_.push_back(xH2O[i]);
    ref_pres_.push_back(pres[i]);
  }
}

MwrAbsorberNH3::MwrAbsorberNH3(MeshBlock *pmb, ParameterInput *pin,
    std::vector<int> imols, Real xHe, Real power):
  Absorber("mw_NH3"), method_(3), xHe_(xHe), power_(power)
{
  imols_ = imols;
  mixrs_ = {1., 1.};

  std::stringstream msg;

  if (imols_.size() != 2) {
    msg << "### FATAL ERROR in MwrAbsorberNH3::MwrAbsorberNH3."
        << std::endl << "Number of dependent molecules is not 2";
    ATHENA_ERROR(msg);
  }

  if ((xHe < 0.) || (xHe > 1.)) {
    msg << "### FATAL ERROR in MwrAbsorberNH3::MwrAbsorberNH3."
        << std::endl << "Value error in molar mixing ratios";
    ATHENA_ERROR(msg);
  }
}

Real MwrAbsorberNH3::GetAttenuation(Real wave1, Real wave2,
    CellVariables const& var) const
{
  Real P = var.w[IPR]/1.E5; // pa -> bar
  Real P_idl = var.w[IPR]/1.E5; // pa -> bar
  Real T = var.w[IDN];
  Real xdry = 1.;
  for (int i = 1; i <= NVAPOR; ++i) xdry -= var.w[i];
  Real XHe = xHe_*xdry;
  Real XH2, XNH3, XH2O;

  if (method_ == 1) {
    XNH3 = var.w[imols_[0]];
    XH2O = xH2O_*xdry;
    XH2 = xdry - XHe - XH2O;
  } else if (method_ == 2) {
    XNH3 = var.w[imols_[0]];
    XH2O = interp1(var.w[IPR], ref_xh2o_.data(), ref_pres_.data(), ref_pres_.size())*xdry;
    XH2 = xdry - XHe - XH2O;
  } else {  // method_ == 3
    XNH3 = var.w[imols_[0]];
    XH2O = var.w[imols_[1]];
    XH2 = xdry - XHe;
  }

  Real abs;

  if (model_name_ == "Bellotti16")
    abs = attenuation_NH3_Bellotti(wave1, P, P_idl, T, XH2, XHe, XNH3, XH2O);
  else if (model_name_ == "BellottiSwitch16")
    abs = attenuation_NH3_Bellotti_switch(wave1, P, P_idl, T, XH2, XHe, XNH3, XH2O);
  else if (model_name_ == "Devaraj")
    abs = attenuation_NH3_Devaraj(wave1, P, P_idl, T, XH2, XHe, XNH3, XH2O);
  else if (model_name_ == "Radtran")
    abs = attenuation_NH3_radtran(wave1, P, T, XH2, XHe, XNH3);
  else // Hanley09
    abs = attenuation_NH3_Hanley(wave1, P, P_idl, T, XH2, XHe, XNH3, XH2O, power_);

  return 100.*abs;  // 1/cm -> 1/m
}
