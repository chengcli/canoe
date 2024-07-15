// cantera
#include <cantera/kinetics.h>
#include <cantera/kinetics/Condensation.h>
#include <cantera/thermo.h>

// canoe
#include <air_parcel.hpp>

// snap
#include "atm_thermodynamics.hpp"

void Thermodynamics::EquilibrateTP(AirParcel* air) const {
  air->ToMassFraction();
  Real pres = air->w[IPR];
  for (int n = IVX; n < IVX + NCLOUD; ++n) air->w[n] = air->c[n - IVX];
  auto& thermo = kinetics_->thermo();
  thermo.setMassFractionsPartial(&air->w[1]);
  thermo.setDensity(air->w[IDN]);
  thermo.setPressure(air->w[IPR]);

  EquilibrateTP();

  for (int n = 0; n < NCLOUD; ++n) air->c[n] = air->w[IVX + n];
  air->w[IPR] = pres;
  air->ToMoleFraction();
}

void Thermodynamics::EquilibrateTP(Tensor const& w) {
  size_t size = w.sizes();

  thermo.setMassFractionsPartial(&w(1, k, j, i), /*stride=*/size);
  thermo.setDensity(w(IDN, k, j, i));
  thermo.setPressure(w(IPR, k, j, i));

  kinetics_->setQuantityMoleFraction();

  std::vector<Real> rates(kinetics_->nTotalSpecies());
  std::vector<Real> mfrac(kinetics_->nTotalSpecies());

  auto& thermo = kinetics_->thermo();
  Real temp = thermo.temperature();
  Real pres = thermo.pressure();

  int iter = 0, max_iter = 3;

  while (iter++ < max_iter) {
    // get mole fraction
    kinetics_->getActivityConcentrations(mfrac.data());
    kinetics_->getNetProductionRates(rates.data());

    // update mole fraction
    for (size_t i = 1; i < rates.size(); ++i) {
      mfrac[i] += rates[i];
    }

    thermo.setMoleFractions(mfrac.data());
    thermo.setTemperature(temp);
    thermo.setPressure(pres);
  }
}
