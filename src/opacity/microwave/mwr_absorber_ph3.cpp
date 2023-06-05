// C/C++
#include <stdexcept>
#include <sstream>

// canoe
#include <configure.hpp>

// athena
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// climath
#include <climath/interpolation.h>

// snap
#include <snap/cell_variables.hpp>

// opacity
#include "mwr_absorbers.hpp"
#include "absorption_functions.hpp"

MwrAbsorberPH3::MwrAbsorberPH3(MeshBlock *pmb, ParameterInput *pin,
    int imol, Real xHe):
  Absorber(pmb, pin, "", "mw_PH3"), xHe_(xHe), method_(1)
{
  imols_ = {imol};
  mixrs_ = {1.};

  std::stringstream msg;
  if ((xHe_ < 0.) || (xHe_ > 1.)) {
    msg << "### FATAL ERROR in MwrAbsorberPH3::MwrAbsorberPH3."
        << std::endl << "Value error in molar mixing ratios";
    throw std::runtime_error(msg.str().c_str());
  }
}

MwrAbsorberPH3::MwrAbsorberPH3(MeshBlock *pmb, ParameterInput *pin,
    Real xHe, Real *xPH3, Real *pres, int np):
  Absorber(pmb, pin, "", "mw_PH3"), xHe_(xHe), method_(2)
{
  std::stringstream msg;
  for (int i = 0; i < np; ++i) {
    if ((xPH3[i] < 0.) || (xPH3[i] > 1.)) {
      msg << "### FATAL ERROR in MwrAbsorberNH3::MwrAbsorberNH3."
          << std::endl << "Value error in molar mixing ratios";
      throw std::runtime_error(msg.str().c_str());
    }
    ref_xph3_.push_back(xPH3[i]);
    ref_pres_.push_back(pres[i]);
  }
}

Real MwrAbsorberPH3::getAttenuation(Real wave1, Real wave2,
    CellVariables const& var) const
{
  Real P = var.w[IPR]/1.E5; // pa -> bar
  Real T = var.w[IDN];
  Real xdry = 1.;
  for (int i = 1; i <= NVAPOR; ++i) xdry -= var.w[i];
  Real XHe = xHe_*xdry;
  Real XH2, XPH3;

  if (method_ == 1) {
    XPH3 = var.w[imols_[0]];
    XH2 = xdry - XHe;
  } else {  // method_ == 2
    XPH3 = interp1(var.w[IPR], ref_xph3_.data(), ref_pres_.data(), ref_pres_.size())*xdry;;
    XH2 = xdry - XHe - XPH3;
  }

  Real abs;

  if (model_name_ == "Radtran")
    abs = absorption_coefficient_PH3_radtran(wave1, P, T, XH2, XHe, XPH3);
  else // Hoffman
    abs = absorption_coefficient_PH3_Hoffman(wave1, P, T, XH2, XHe, XPH3);

  return 100.*abs;  // 1/cm -> 1/m
}
