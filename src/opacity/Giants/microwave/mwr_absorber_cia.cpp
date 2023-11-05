// canoe
#include <air_parcel.hpp>
#include <configure.hpp>

// application
#include <application/exceptions.hpp>

// opacity
#include "absorption_functions.hpp"
#include "mwr_absorbers.hpp"

namespace GiantPlanets {

void MwrAbsorberCIA::CheckFail() const {
  if (HasPar("xHe")) {
    throw NotFoundError("MwrAbsorberCIA", "parameter 'xHe'");
  }

  if (HasPar("xCH4")) {
    throw NotFoundError("MwrAbsorberCIA", "parameter 'xCH4'");
  }

  if (HasPar("mix")) {
    throw NotFoundError("MwrAbsorberCIA", "parameter 'mix'");
  }
}

Real MwrAbsorberCIA::GetAttenuation(Real wave1, Real wave2,
                                    AirParcel const& var) const {
  Real P = var.w[IPR] / 1.E5;  // pa -> bar
  Real T = var.w[IDN];
  Real xdry = 1.;
  for (int i = 1; i <= NVAPOR; ++i) xdry -= var.w[i];
  Real XHe = GetPar<Real>("xHe") * xdry;
  Real XCH4 = GetPar<Real>("xCH4") * xdry;
  Real XH2 = (1. - GetPar<Real>("xHe") - GetPar<Real>("xCH4")) * xdry;
  Real mix = GetPar<Real>("mix");
  Real wave = (wave1 + wave2) / 2.;

  // 1/cm -> 1/m
  return 100. * attenuation_CIA(wave, P, T, XH2, XHe, XCH4, mix);
}

}  // namespace GiantPlanets
