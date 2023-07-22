// canoe
#include <configure.hpp>
#include <variable.hpp>

// opacity
#include "absorption_functions.hpp"
#include "mwr_absorbers.hpp"

namespace GiantPlanets {

Real MwrAbsorberH2S::GetAttenuation(Real wave1, Real wave2,
                                    Variable const& var) const {
  // adapted by cli (Cheng Li), Aug 30
  Real P = var.w[IPR] / 1.E5;
  Real T = var.w[IDN];
  Real xdry = 1.;
  for (int i = 1; i <= NVAPOR; ++i) xdry -= var.w[i];

  Real XHe = params_.at("xHe") * xdry;
  Real XH2S = var.w[imols_[0]];
  Real XH2 = xdry - XHe;
  Real wave = (wave1 + wave2) / 2.;

  // 1/cm -> 1/m
  return 100. * attenuation_H2S_Hofstadter(wave, P, T, XH2, XHe, XH2S);
}

}  // namespace GiantPlanets
