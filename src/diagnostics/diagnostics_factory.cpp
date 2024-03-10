// C/C++
#include <cstring>

// athena
#include <athena/parameter_input.hpp>

// diagnostics
#include "diagnostics.hpp"

DiagnosticsContainer DiagnosticsFactory::CreateFrom(MeshBlock *pmb,
                                                    ParameterInput *pin) {
  DiagnosticsContainer diag;

  char cstr[80];
  std::string diag_names = pin->GetOrAddString("problem", "diagnostics", "");
  std::strcpy(cstr, diag_names.c_str());
  char *p = std::strtok(cstr, " ,");

  while (p != NULL) {
    std::string name(p);
    if (name == "div") {  // 1.
      diag.push_back(std::make_shared<Divergence>(pmb));
    } else if (name == "curl") {  // 2.
      diag.push_back(std::make_shared<Curl>(pmb));
    } else if (name == "b") {  // 3.
      diag.push_back(std::make_shared<Buoyancy>(pmb));
    } else if (name == "mean") {  // 4.
      diag.push_back(std::make_shared<HydroMean>(pmb));
    } else if (name == "div_h") {  // 5.
      diag.push_back(std::make_shared<HorizontalDivergence>(pmb));
    } else if (name == "tempa") {  // 6.
      diag.push_back(std::make_shared<TemperatureAnomaly>(pmb));
    } else if (name == "presa") {  // 7.
      diag.push_back(std::make_shared<PressureAnomaly>(pmb));
    } else if (name == "radflux") {  // 8.
      diag.push_back(std::make_shared<RadiativeFlux>(pmb));
    } else if (name == "hydroflux") {  // 9.
      diag.push_back(std::make_shared<HydroFlux>(pmb));
      /*} else if (name == "eddyflux") { // 6.
        diag.push_back(std::make_shared<EddyFlux>(pmb));
      } else if (name == "am") { // 11.
        diag.push_back(std::make_shared<AngularMomentum>(pmb));
      } else if (name == "eke") { // 12.
        diag.push_back(std::make_shared<EddyKineticEnergy>(pmb));
      } else if (name == "tendency") { // 13.
        diag.push_back(std::make_shared<Tendency>(pmb)); */
    } else {
      throw std::invalid_argument("Invalid diagnostic variable: " + name);
    }
    p = std::strtok(NULL, " ,");
  }

  return diag;
}
