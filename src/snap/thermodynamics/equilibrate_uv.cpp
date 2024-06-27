// cantera
#include <cantera/kinetics.h>
#include <cantera/kinetics/Condensation.h>
#include <cantera/thermo.h>

// canoe
#include <air_parcel.hpp>

// snap
#include "atm_thermodynamics.hpp"

void Thermodynamics::EquilibrateUV(AirParcel* air) const {
  air->ToMassFraction();
  Real pres = air->w[IPR];
  for (int n = IVX; n < IVX + NCLOUD; ++n) air->w[n] = air->c[n - IVX];
  auto& thermo = kinetics_->thermo();
  thermo.setMassFractionsPartial(&air->w[1]);
  thermo.setDensity(air->w[IDN]);
  thermo.setPressure(air->w[IPR]);

  EquilibrateUV();

  for (int n = 0; n < NCLOUD; ++n) air->c[n] = air->w[IVX + n];
  air->w[IPR] = pres;
  air->ToMoleFraction();
}

void Thermodynamics::EquilibrateUV() const {
  auto kin = std::static_pointer_cast<Cantera::Condensation>(kinetics_);
  kin->setQuantityConcentration();

  Eigen::VectorXd rates(kin->nTotalSpecies());
  Eigen::VectorXd conc(kin->nTotalSpecies());
  Eigen::VectorXd intEng(kin->nTotalSpecies());
  Eigen::VectorXd cv(kin->nTotalSpecies());

  auto& thermo = kin->thermo();

  int iter = 0, max_iter = 3;
  while (iter++ < max_iter) {
    /*std::cout << "#############" << std::endl;
    std::cout << "Iteration " << iter << std::endl;*/

    Real temp = thermo.temperature();
    // Real pres = thermo.pressure();

    // get concentration
    kin->getActivityConcentrations(conc.data());
    thermo.getIntEnergy_RT(intEng.data());
    thermo.getCv_R(cv.data());

    Real cc = conc.dot(cv);
    Real uc = conc.dot(intEng);

    /*std::cout << "internal energy 1 = " << uc * temp << std::endl;

    // print initial conditions
    std::cout << "Concentration" << std::endl;
    for (size_t i = 0; i < conc.size(); ++i) {
      std::cout << "Concentration" << i << ": " << conc[i] << std::endl;
    }*/

    kin->getNetProductionRates(rates.data());
    // std::cout << "rates = " << rates << std::endl;

    // update concentrations
    for (size_t i = 1; i < rates.size(); ++i) {
      conc(i) += rates(i);
    }

    // update temperature
    auto dT = -temp * rates.dot(intEng) / conc.dot(cv);
    // std::cout << "dT = " << dT << std::endl;

    /* print initial conditions
    std::cout << "Concentration" << std::endl;
    for (size_t i = 0; i < conc.size(); ++i) {
      std::cout << "Concentration" << i << ": " << conc[i] << std::endl;
    }*/

    thermo.setConcentrationsNoNorm(conc.data());
    thermo.setTemperature(temp + dT);

    /*thermo.getIntEnergy_RT(intEng.data());
    thermo.getCv_R(cv.data());
    std::cout << "internal energy 2 = "
              << conc.dot(intEng) * thermo.temperature() << std::endl;

    std::cout << "T = " << thermo.temperature() << std::endl;
    std::cout << "P = " << thermo.pressure() << std::endl;
    std::cout << "D = " << thermo.density() << std::endl;*/
  }
}
