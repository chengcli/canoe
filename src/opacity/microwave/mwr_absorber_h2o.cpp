// C/C++
#include <sstream>
#include <stdexcept>

// canoe
#include <configure.hpp>

// athena
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// snap
#include <snap/cell_variables.hpp>

// opacity
#include "mwr_absorbers.hpp"
#include "absorption_functions.hpp"

MwrAbsorberH2O::MwrAbsorberH2O(MeshBlock *pmb, ParameterInput *pin,
    int imol, Real xHe, Real scale):
  Absorber("mw_H2O"), xHe_(xHe), scale_(scale)
{
  imols_ = {imol};
  mixrs_ = {1.};
  std::stringstream msg;

  if ((xHe_ < 0.) || (xHe_ > 1.)) {
    msg << "### FATAL ERROR in MwrAbsorberPH3::MwrAbsorberH2O."
        << std::endl << "Value error in molar mixing ratios";
    ATHENA_ERROR(msg);
  }
}

Real MwrAbsorberH2O::GetAttenuation(Real wave1, Real wave2,
    CellVariables const& var) const
{
  Real P = var.w[IPR]/1.E5; // pa -> bar
  Real T = var.w[IDN];
  Real xdry = 1.;
  for (int i = 1; i <= NVAPOR; ++i) xdry -= var.w[i];
  Real XHe = xHe_*xdry;
  Real XH2 = xdry - XHe;
  Real XH2O = var.w[imols_[0]];

  Real abs;

  if (model_name_ == "deBoer")
    abs = attenuation_H2O_deBoer(wave1, P, T, XH2, XHe, XH2O);
  else if (model_name_ == "Waters")
    abs = attenuation_H2O_Waters(wave1, P, T, XH2, XHe, XH2O);
  else if (model_name_ == "Goodman")
    abs = attenuation_H2O_Goodman(wave1, P, T, XH2, XHe, XH2O);
  else // Karpowicz
    abs = attenuation_H2O_Karpowicz(wave1, P, T, XH2, XHe, XH2O, scale_);

  return 100.*abs;  // 1/cm -> 1/m
}
