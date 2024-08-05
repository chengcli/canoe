// Eigen
#include <Eigen/Core>

// cantera
#include <cantera/kinetics.h>
#include <cantera/kinetics/Condensation.h>
#include <cantera/thermo.h>

// snap
#include "thermodynamics.hpp"

void Thermodynamics::EquilibrateUV() const {
  kinetics_->setQuantityConcentration();

  Eigen::VectorXd rates(Size);
  Eigen::VectorXd conc(Size);
  Eigen::VectorXd intEng(Size);
  Eigen::VectorXd cv(Size);

  auto& thermo = kinetics_->thermo();

  for (int iter = 0; iter < 3; ++iter) {
    // std::cout << "#############" << std::endl;
    // std::cout << "Iteration " << iter << std::endl;

    kinetics_->getNetProductionRates(rates.data());

    // get concentration
    kinetics_->getActivityConcentrations(conc.data());
    thermo.getIntEnergy_RT(intEng.data());
    thermo.getCv_R(cv.data());

    Real cc = conc.dot(cv);
    Real uc = conc.dot(intEng);

    /* print initial conditions
    std::cout << "internal energy 1 = " << uc * temp << std::endl;
    std::cout << "Concentration" << std::endl;
    for (size_t i = 0; i < conc.size(); ++i) {
      std::cout << "Concentration" << i << ": " << conc[i] << std::endl;
    }*/

    // update concentrations
    for (size_t i = 1; i < rates.size(); ++i) {
      conc(i) += rates(i);
    }

    // update temperature
    Real temp = thermo.temperature();
    auto dT = -temp * rates.dot(intEng) / conc.dot(cv);

    /* print initial conditions
    std::cout << "rates = " << rates << std::endl;
    std::cout << "dT = " << dT << std::endl;
    std::cout << "Concentration" << std::endl;
    for (size_t i = 0; i < conc.size(); ++i) {
      std::cout << "Concentration" << i << ": " << conc[i] << std::endl;
    }*/

    thermo.setConcentrationsNoNorm(conc.data());
    thermo.setTemperature(temp + dT);

    if (std::abs(dT) < 0.1) break;

    /*thermo.getIntEnergy_RT(intEng.data());
    thermo.getCv_R(cv.data());
    std::cout << "internal energy 2 = "
              << conc.dot(intEng) * thermo.temperature() << std::endl;
    std::cout << "T = " << thermo.temperature() << std::endl;
    std::cout << "P = " << thermo.pressure() << std::endl;
    std::cout << "D = " << thermo.density() << std::endl;*/
  }
}
