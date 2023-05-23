/** @file ionization.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Saturday Dec 04, 2021 23:30:27 EST
 * @bug No known bugs.
 */

// Athena++ header
#include "../radiation/radiation.hpp"
#include "thermodynamics.hpp"

Real saha_ionization_electron_density(Real T, Real num, Real ion_ev) {
  num *= 1.E-6;  // m^{-3} -> cm^{-3}
  Real erg2ev = 6.242E11;
  Real kBoltz_ev = Thermodynamics::kBoltz_cgs * erg2ev;
  // Real kBoltz_ev = 8.617E-5;  // ev/K
  Real hPlanck_ev = Radiation::hPlanck_cgs * erg2ev;  // ev.s
  // Real hPlanck_ev = 4.136E-15; // ev.s
  Real me_cgs = 9.11E-28;    // g
  Real kBT = kBoltz_ev * T;  // ev
  // Real kBT = T/300.*0.0256; // ev
  Real hbar = hPlanck_ev / (2. * M_PI);
  // Real hbarc = 1973E-8; // ev.cm
  Real hbarc = hbar * Radiation::cLight_cgs;  // ev.cm
  Real mc2 = 0.511E6;                         // ev
  Real debroglie_thermal = sqrt(2. * M_PI * hbarc * hbarc / (mc2 * kBT));  // cm
  Real coeff = 2. / pow(debroglie_thermal, 3) * exp(-ion_ev / kBT);  // cm^{-3}
  return 1.E6 * (-coeff + sqrt(coeff * coeff + 4. * coeff * num)) /
         2.;  // cm^{-3} -> m^{-3}
}
