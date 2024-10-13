// cantera
#include <cantera/kinetics.h>
#include <cantera/kinetics/Condensation.h>
#include <cantera/thermo.h>

// snap
#include "thermodynamics.hpp"

void Thermodynamics::EquilibrateTP(Real temp, Real pres) const {
  kinetics_->setQuantityMoleFraction();

  std::array<Real, Size> rates;
  std::array<Real, Size> xfrac;

  auto& thermo = kinetics_->thermo();

  thermo.setTemperature(temp);
  thermo.setPressure(pres);

  int iter = 0, max_iter = 3;
  while (iter++ < max_iter) {
    // get mole fraction
    kinetics_->getNetProductionRates(rates.data());

    Real max_abs_rate = 0.;
    for (size_t i = 1; i < Size; ++i) {
      if (std::abs(rates[i]) > max_abs_rate) {
        max_abs_rate = std::abs(rates[i]);
      }
    }
    if (max_abs_rate < 1.E-8) break;

    thermo.getMoleFractions(xfrac.data());

    /*std::cout << "iter: " << iter << std::endl;
    for (size_t i = 0; i < Size; ++i) {
      std::cout << xfrac[i] << ", " << rates[i] << ", ";
    }
    std::cout << std::endl;*/

    // update mole fraction
    for (size_t i = 1; i < Size; ++i) {
      xfrac[i] += rates[i];
    }

    thermo.setMoleFractions(xfrac.data());

    // updates density
    thermo.setTemperature(temp);
    thermo.setPressure(pres);
  }
}
