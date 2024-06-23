// cantera
#include <cantera/kinetics.h>
#include <cantera/kinetics/Condensation.h>
#include <cantera/thermo.h>

// canoe
#include <air_parcel.hpp>

// snap
#include "atm_thermodynamics.hpp"

void Thermodynamics::EquilibrateTP(AirParcel* air) const {}

void Thermodynamics::EquilibrateTP() const {
  std::static_pointer_cast<Cantera::Condensation>(kinetics_)
      ->setQuantityMoleFraction();

  std::vector<Real> rates(kinetics_->nTotalSpecies());
  std::vector<Real> mfrac(kinetics_->nTotalSpecies());

  auto& thermo = kinetics_->thermo();
  Real temp = thermo.temperature();
  Real pres = thermo.pressure();

  int iter = 0, max_iter = 5;

  while (iter++ < max_iter) {
    std::cout << "Iteration " << iter << std::endl;

    // get mole fraction
    kinetics_->getActivityConcentrations(mfrac.data());

    // print initial conditions
    std::cout << "Mole fractions" << std::endl;
    for (size_t i = 0; i < mfrac.size(); ++i) {
      std::cout << "Mole fraction " << i << ": " << mfrac[i] << std::endl;
    }

    kinetics_->getNetProductionRates(rates.data());
    for (size_t i = 0; i < rates.size(); ++i) {
      std::cout << "NET Production " << i << ": " << rates[i] << std::endl;
    }

    for (size_t i = 1; i < rates.size(); ++i) {
      mfrac[i] += rates[i];
    }

    thermo.setMoleFractions(mfrac.data());
    thermo.setTemperature(temp);
    thermo.setPressure(pres);

    // print again
    std::cout << "After Mole fractions" << std::endl;
    for (size_t i = 0; i < mfrac.size(); ++i) {
      std::cout << "Mole fraction " << i << ": " << mfrac[i] << std::endl;
    }

    std::cout << "T = " << thermo.temperature() << std::endl;
    std::cout << "P = " << thermo.pressure() << std::endl;
    std::cout << "D = " << thermo.density() << std::endl;
  }
}
