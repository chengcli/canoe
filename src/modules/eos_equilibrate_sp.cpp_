// cantera
#include <cantera/kinetics.h>
#include <cantera/kinetics/Condensation.h>
#include <cantera/thermo.h>

// application
#include <application/exceptions.hpp>

// climath
#include <climath/root.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

void Thermodynamics::EquilibrateSP(double P) const {
  std::static_pointer_cast<Cantera::Condensation>(kinetics_)
      ->setQuantityMoleFraction();

  auto& thermo = kinetics_->thermo();

  double entropy = thermo.entropy_mole();
  Real temp = thermo.temperature();

  int err = root(temp / 2., temp * 2., 1.E-8, &temp, [&](Real T) {
    thermo.setTemperature(T);
    thermo.setPressure(P);
    EquilibrateTP();
    return thermo.entropy_mole() - entropy;
  });
  if (err)
    throw RuntimeError("Thermodynamics::EquilibrateSP",
                       "EquilibrateSP doesn't converge");

  thermo.setTemperature(temp);
  thermo.setPressure(P);
  EquilibrateTP();
}
