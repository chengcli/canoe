// canoe
#include <air_parcel.hpp>
#include <configure.hpp>

// opacity
#include "absorption_functions.hpp"
#include "mwr_absorbers.hpp"

namespace GiantPlanets {

Real MwrAbsorberElectron::GetAttenuation(Real wave1, Real wave2,
                                         AirParcel const& qfrac) const {
  Real P = qfrac.w[IPR] / 1.E5;  // pa -> bar
  Real T = qfrac.w[IDN];

  Real abs;
  Real wave = (wave1 + wave2) / 2.;

  AirParcel var(qfrac);
  var.ToMoleConcentration();

  if (model_name_ == "Reference") {
    abs = attenuation_freefree_Reference(wave, P, T);
  } else if (model_name_ == "ChengLi") {
    abs = attenuation_freefree_Chengli(wave, P, T);
  } else {  // AppletonHartree
    abs = attenuation_appleton_hartree_nomag(wave, P, T,
                                             var.w[GetSpeciesIndex(0)]);
  }

  return 100. * abs;  // 1/cm -> 1/m
}

}  // namespace GiantPlanets
