// cantera
#include <cantera/kinetics.h>
#include <cantera/kinetics/Condensation.h>
#include <cantera/thermo.h>

// snap
#include "thermodynamics.hpp"

void Thermodynamics::EquilibrateTP() const {
  kinetics_->setQuantityMoleFraction();

  std::vector<Real> rates(Size);
  std::vector<Real> xfrac(Size);

  auto& thermo = kinetics_->thermo();

  Real temp = thermo.temperature();
  Real pres = thermo.pressure();

  int iter = 0, max_iter = 3;

  while (iter++ < max_iter) {
    // get mole fraction
    kinetics_->getActivityConcentrations(xfrac.data());
    kinetics_->getNetProductionRates(rates.data());
    /*for (size_t i = 0; i < rates.size(); ++i) {
      std::cout << rates[i] << ", ";
    }
    std::cout << std::endl;*/

    // update mole fraction
    for (size_t i = 1; i < rates.size(); ++i) {
      xfrac[i] += rates[i];
    }

    thermo.setMoleFractions(xfrac.data());
    thermo.setTemperature(temp);
    thermo.setPressure(pres);
  }
}
