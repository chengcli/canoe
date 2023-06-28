// C/C++
#include <sstream>
#include <stdexcept>

// canoe
#include <configure.hpp>
#include <variable.hpp>

// opacity
#include "absorption_functions.hpp"
#include "mwr_absorbers.hpp"

namespace GiantPlanets {

Real MwrAbsorberCIA::GetAttenuation(Real wave1, Real wave2,
                                    Variable const& var) const {
  Real P = var.w[IPR] / 1.E5;  // pa -> bar
  Real T = var.w[IDN];
  Real xdry = 1.;
  for (int i = 1; i <= NVAPOR; ++i) xdry -= var.w[i];
  Real XHe = params_.at("xHe") * xdry;
  Real XCH4 = params_.at("xCH4") * xdry;
  Real XH2 = (1. - params_.at("xHe") - params_.at("xCH4")) * xdry;
  Real fequal = params_.at("fequal");

  // 1/cm -> 1/m
  return 100. * absorption_coefficient_CIA(wave1, P, T, XH2, XHe, XCH4, fequal);
}

}  // namespace GiantPlanets
