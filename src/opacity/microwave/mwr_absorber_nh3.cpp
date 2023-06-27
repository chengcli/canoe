// C/C++
#include <sstream>
#include <stdexcept>

// canoe
#include <configure.hpp>

// athena
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// climath
#include <climath/interpolation.h>

// application
#include <application/exceptions.hpp>

// snap
#include <snap/cell_variables.hpp>

// opacity
#include "mwr_absorbers.hpp"
#include "absorption_functions.hpp"


Real MwrAbsorberNH3::GetAttenuation(Real wave1, Real wave2,
    CellVariables const& var) const
{
  Real P = var.w[IPR]/1.E5; // pa -> bar
  Real P_idl = var.w[IPR]/1.E5; // pa -> bar
  Real T = var.w[IDN];
  Real xdry = 1.;
  for (int i = 1; i <= NVAPOR; ++i) xdry -= var.w[i];

  Real XHe = params_.at("xHe")*xdry;
  Real XH2 = xdry - XHe;
  Real XNH3 = var.w[imols_[0]];
  Real XH2O = var.w[imols_[1]];

  Real abs;

  if (model_name_ == "Bellotti16") {
    abs = attenuation_NH3_Bellotti(wave1, P, P_idl, T, XH2, XHe, XNH3, XH2O);
  } else if (model_name_ == "BellottiSwitch16") {
    abs = attenuation_NH3_Bellotti_switch(wave1, P, P_idl, T, XH2, XHe, XNH3, XH2O);
  } else if (model_name_ == "Devaraj") {
    abs = attenuation_NH3_Devaraj(wave1, P, P_idl, T, XH2, XHe, XNH3, XH2O);
  } else if (model_name_ == "Radtran") {
    abs = attenuation_NH3_radtran(wave1, P, T, XH2, XHe, XNH3);
  } else if (model_name_ == "Hanley09") {
    abs = attenuation_NH3_Hanley(wave1, P, P_idl, T, XH2, XHe, XNH3, XH2O, params_.at("power"));
  } else {
    throw NotFoundError("MwrAbsorberNH3::GetAttenuation: " + model_name_ );
  }

  return 100.*abs;  // 1/cm -> 1/m
}
