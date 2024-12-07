// Eigen
#include <Eigen/Core>

// cantera
#include <cantera/kinetics.h>
#include <cantera/kinetics/Condensation.h>
#include <cantera/thermo.h>

// snap
#include "thermodynamics.hpp"

void Thermodynamics::EquilibrateUV(Real dt) const {
  kinetics_->setQuantityConcentration(dt);

  Eigen::VectorXd rates(Size);
  Eigen::VectorXd conc(Size);
  Eigen::VectorXd conc0(Size);
  Eigen::VectorXd intEng(Size);
  Eigen::VectorXd cv(Size);

  auto& thermo = kinetics_->thermo();

  for (int iter = 0; iter < 10; ++iter) {
    // std::cout << "#############" << std::endl;
    // std::cout << "Iteration " << iter << std::endl;

    kinetics_->getNetProductionRates(rates.data());

    // get concentration
    thermo.getConcentrations(conc.data());
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
    conc0 = conc;
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

    if (temp + dT < 0) {
      for (size_t i = 1; i < rates.size(); ++i) {
        rates(i) /= 2.;
        conc(i) = conc0(i) + rates(i);
      }
      dT = -temp * rates.dot(intEng) / conc.dot(cv);
    }

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
