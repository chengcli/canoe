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

  int iter = 0, max_iter = 3;
  while (iter++ < max_iter) {
    /*std::cout << "#############" << std::endl;
    std::cout << "Iteration " << iter << std::endl;*/

    Real temp = thermo.temperature();
    Real pres = thermo.pressure();

    kinetics_->getNetProductionRates(rates.data());

    /*std::cout << "temp = " << temp << std::endl;
    std::cout << "pres = " << pres << std::endl;
    std::cout << "rates = " << rates << std::endl;*/

    Real max_abs_rate = 0.;
    for (size_t i = 1; i < Size; ++i) {
      if (std::abs(rates[i]) > max_abs_rate) {
        max_abs_rate = std::abs(rates[i]);
      }
    }
    if (max_abs_rate < 1.E-8) break;

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
    auto dT = -temp * rates.dot(intEng) / conc.dot(cv);

    /* print initial conditions
    std::cout << "dT = " << dT << std::endl;
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
