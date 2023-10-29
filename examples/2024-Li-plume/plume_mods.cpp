// application
#include <application/application.hpp>

// canoe
#include <index_map.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

double sat_vapor_p_H2O_liquid_Ideal(double T) {
  double betal = 24.845, gammal = 4.986009, tr = 273.16, pr = 611.7;
  return SatVaporPresIdeal(T / tr, pr, betal, gammal);
}

void Thermodynamics::enrollVaporFunctions() { enrollVaporFunctionH2O(); }

// water svp
void Thermodynamics::enrollVaporFunctionH2O() {
  Application::Logger app("snap");
  app->Log("Enrolling H2O vapor pressures");

  auto pindex = IndexMap::GetInstance();
  int iH2O = pindex->GetVaporId("H2O");

  svp_func1_[iH2O][0] = [](AirParcel const &qfrac, int, int) {
    return sat_vapor_p_H2O_liquid_Ideal(qfrac.w[IDN]);
  };
}
