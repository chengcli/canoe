// application
#include <application/application.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

/*double sat_vapor_p_H2O(double T) {
  double betal = 24.845, gammal = -2.1735, tr = 273.16, pr = 611.7;
  return SatVaporPresIdeal(T / tr, pr, betal, gammal);
}

// water svp
void enroll_vapor_Function_H2O(Thermodynamics::SVPFunc1Container &svp_func1) {
  Application::Logger app("snap");
  app->Log("Enrolling H2O vapor pressures");

  auto pindex = IndexMap::GetInstance();
  int iH2O = pindex->GetVaporId("H2O");

  svp_func1[iH2O][0] = [](AirParcel const &qfrac, int, int) {
    return sat_vapor_p_H2O(qfrac.w[IDN]);
  };
}

void Thermodynamics::enrollVaporFunctions() {
  enroll_vapor_Function_H2O(svp_func1_);
}*/
