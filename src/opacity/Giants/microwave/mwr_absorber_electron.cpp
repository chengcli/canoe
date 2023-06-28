// canoe
#include <configure.hpp>
#include <variable.hpp>

// opacity
#include "absorption_functions.hpp"
#include "mwr_absorbers.hpp"

namespace GiantPlanets {

Real MwrAbsorberElectron::GetAttenuation(Real wave1, Real wave2,
                                         Variable const& var) const {
  Real P = var.w[IPR] / 1.E5;  // pa -> bar
  Real T = var.w[IDN];

  Real abs;
  Real wave = (wave1 + wave2) / 2.;

  if (model_name_ == "Reference") {
    abs = attenuation_freefree_Reference(wave, P, T);
  } else if (model_name_ == "ChengLi") {
    abs = attenuation_freefree_Chengli(wave, P, T);
  } else {  // AppletonHartree
    abs = attenuation_appleton_hartree_nomag(wave, P, T, var.s[imols_[0]]);
  }

  return 100. * abs;  // 1/cm -> 1/m
}

}  // namespace GiantPlanets
