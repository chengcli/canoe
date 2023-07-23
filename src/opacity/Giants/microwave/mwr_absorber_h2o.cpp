// canoe
#include <configure.hpp>
#include <variable.hpp>

// application
#include <application/exceptions.hpp>

// opacity
#include "absorption_functions.hpp"
#include "mwr_absorbers.hpp"

namespace GiantPlanets {

MwrAbsorberH2O::MwrAbsorberH2O(SpeciesNames const& species, ParameterMap params)
    : Absorber("radio-H2O", species, params) {
  if (!params_.count("xHe")) {
    throw NotFoundError("MwrAbsorberH2O", "parameter 'xHe'");
  }

  if (!params_.count("scale")) {
    throw NotFoundError("MwrAbsorberH2O", "parameter 'scale'");
  }
}

Real MwrAbsorberH2O::GetAttenuation(Real wave1, Real wave2,
                                    Variable const& var) const {
  Real P = var.w[IPR] / 1.E5;  // pa -> bar
  Real T = var.w[IDN];
  Real xdry = 1.;
  for (int i = 1; i <= NVAPOR; ++i) xdry -= var.w[i];
  Real XHe = params_.at("xHe") * xdry;
  Real XH2 = xdry - XHe;
  Real XH2O = var.w[imols_[0]];

  Real abs;
  Real wave = (wave1 + wave2) / 2.;

  if (model_name_ == "deBoer") {
    abs = attenuation_H2O_deBoer(wave, P, T, XH2, XHe, XH2O);
  } else if (model_name_ == "Waters") {
    abs = attenuation_H2O_Waters(wave, P, T, XH2, XHe, XH2O);
  } else if (model_name_ == "Goodman") {
    abs = attenuation_H2O_Goodman(wave, P, T, XH2, XHe, XH2O);
  } else {  // Karpowicz
    abs = attenuation_H2O_Karpowicz(wave, P, T, XH2, XHe, XH2O,
                                    params_.at("scale"));
  }

  return 100. * abs;  // 1/cm -> 1/m
}

}  // namespace GiantPlanets
