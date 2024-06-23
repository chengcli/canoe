// cantera
#include <cantera/kinetics.h>
#include <cantera/thermo.h>

// canoe
#include <air_parcel.hpp>

// snap
#include "atm_thermodynamics.hpp"

void Thermodynamics::EquilibrateTP(AirParcel *air) const {}

void Thermodynamics::EquilibrateTP() const {
  std::vector<Real> rates(kinetics_->nTotalSpecies());
  std::vector<Real> conc(kinetics_->nTotalSpecies());

  // save initial concentrations
  for (size_t k = 0; k < kinetics_->nPhases(); ++k) {
    auto &thermo = kinetics_->thermo(k);
    size_t n = kinetics_->kineticsSpeciesIndex(0, k);
    thermo.getConcentrations(&conc[n]);
  }
  Real pres0 = kinetics_->thermo(1).pressure();

  int iter = 0;
  int max_iter = 10;
  while (iter++ < max_iter) {
    std::cout << "Iteration " << iter << std::endl;

    kinetics_->getNetProductionRates(rates.data());
    for (size_t i = 0; i < rates.size(); ++i) {
      std::cout << "NET Production " << i << ": " << rates[i] << std::endl;
    }

    // update concentrations
    for (size_t i = 1; i < rates.size(); ++i) {
      conc[i] += rates[i];
    }

    for (size_t k = 0; k < kinetics_->nPhases(); ++k) {
      auto &thermo = kinetics_->thermo(k);
      size_t n = kinetics_->kineticsSpeciesIndex(0, k);
      thermo.setConcentrationsNoNorm(&conc[n]);
      // thermo.setPressure(49657.6);
      // thermo.setPressure(1.e5);
    }

    // volumne change
    Real x = pres0 / kinetics_->thermo(1).pressure();
    for (size_t i = 0; i < conc.size(); ++i) {
      conc[i] *= x;
    }

    for (size_t k = 0; k < kinetics_->nPhases(); ++k) {
      auto &thermo = kinetics_->thermo(k);
      size_t n = kinetics_->kineticsSpeciesIndex(0, k);
      thermo.setConcentrationsNoNorm(&conc[n]);
      // thermo.setPressure(49657.6);
      // thermo.setPressure(1.e5);
    }

    // print again
    std::cout << "Concentrations" << std::endl;
    for (size_t i = 0; i < conc.size(); ++i) {
      std::cout << "Concentration " << i << ": " << conc[i] << std::endl;
    }
    std::cout << "P = " << kinetics_->thermo(1).pressure() << std::endl;
  }
}
