// canoe
#include <configure.hpp>

// climath
#include <climath/interpolation.h>

// snap
#include <snap/variable.hpp>

// opacity
#include "absorption_functions.hpp"
#include "mwr_absorbers.hpp"

namespace GiantPlanets {

Real MwrAbsorberPH3::GetAttenuation(Real wave1, Real wave2,
                                    Variable const& var) const {
  Real P = var.w[IPR] / 1.E5;  // pa -> bar
  Real T = var.w[IDN];
  Real xdry = 1.;
  for (int i = 1; i <= NVAPOR; ++i) xdry -= var.w[i];

  Real XHe = params_.at("xHe") * xdry;
  Real XPH3 = var.w[imols_[0]];
  Real XH2 = xdry - XHe;

  Real abs;

  if (model_name_ == "Radtran") {
    abs = absorption_coefficient_PH3_radtran(wave1, P, T, XH2, XHe, XPH3);
  } else {  // Hoffman
    abs = absorption_coefficient_PH3_Hoffman(wave1, P, T, XH2, XHe, XPH3);
  }

  return 100. * abs;  // 1/cm -> 1/m
}

}  // namespace GiantPlanets
