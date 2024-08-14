// application
#include <application/application.hpp>

// canoe
#include <index_map.hpp>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

// water svp
double sat_vapor_p_H2O(double T) {
  double betal = 24.845, gammal = -2.1735, tr = 273.16, pr = 611.7;
  return SatVaporPresIdeal(T / tr, pr, betal, gammal);
}

void enroll_vapor_Function_H2O(Thermodynamics::SVPFunc1Container &svp_func1,
                               std::vector<IndexSet> const &cloud_index_set) {
  Application::Logger app("snap");
  app->Log("Enrolling H2O vapor pressures");

  auto pindex = IndexMap::GetInstance();
  if (!pindex->HasVapor("H2O")) return;
  int iH2O = pindex->GetVaporId("H2O");

  svp_func1[iH2O][0] = [](AirParcel const &qfrac, int, int) {
    return sat_vapor_p_H2O(qfrac.w[IDN]);
  };

  for (int n = 1; n < cloud_index_set[iH2O].size(); ++n) {
    svp_func1[iH2O][n] = NullSatVaporPres1;
  }
}

// ammonia svp
double sat_vapor_p_NH3(double T) {
  double betal = 20.08, gammal = 5.62, betas = 20.64, gammas = 1.43, tr = 195.4,
         pr = 6060.;
  return SatVaporPresIdeal(T / tr, pr, betal, gammal);
}

void enroll_vapor_Function_NH3(Thermodynamics::SVPFunc1Container &svp_func1,
                               std::vector<IndexSet> const &cloud_index_set) {
  Application::Logger app("snap");
  app->Log("Enrolling NH3 vapor pressures");

  auto pindex = IndexMap::GetInstance();
  if (!pindex->HasVapor("NH3")) return;
  int iNH3 = pindex->GetVaporId("NH3");

  svp_func1[iNH3][0] = [](AirParcel const &qfrac, int, int) {
    return sat_vapor_p_NH3(qfrac.w[IDN]);
  };

  for (int n = 1; n < cloud_index_set[iNH3].size(); ++n) {
    svp_func1[iNH3][n] = NullSatVaporPres1;
  }
}

// ammonium-sulfide svp
double sat_vapor_p_NH4SH(double T) {
  // xiz
  //   double betal = 70.14, gammal = 4.6, tr = 250.000, pr = 36.94;
  //   return SatVaporPresIdeal(T / tr, pr, betal, gammal);

  // lewis
  return pow(10., 14.82 - 4705. / T) * 101325. * 101325.;

  // Umich
  //   double const GOLB2 = (14.83 - (4715.0 / T));
  //   return (pow(10.0, GOLB2)) * 1013250.0 * 1013250.0;
}

void enroll_vapor_Function_NH4SH(Thermodynamics::SVPFunc1Container &svp_func1,
                                 std::vector<IndexSet> const &cloud_index_set) {
  Application::Logger app("snap");
  app->Log("Enrolling NH4SH vapor pressures");

  auto pindex = IndexMap::GetInstance();
  if (!pindex->HasVapor("NH4SH")) return;
  int iNH4SH = pindex->GetVaporId("NH4SH");

  svp_func1[iNH4SH][0] = [](AirParcel const &qfrac, int, int) {
    return sat_vapor_p_NH4SH(qfrac.w[IDN]);
  };

  for (int n = 1; n < cloud_index_set[iNH4SH].size(); ++n) {
    svp_func1[iNH4SH][n] = NullSatVaporPres1;
  }
}

// silicate svp
double sat_vapor_p_mgsio3(double T) {
  double betal = 36.909, gammal = 0., tr = 1700.000, pr = 6.318;
  return SatVaporPresIdeal(T / tr, pr, betal, gammal);
}

void enroll_vapor_Function_mgsio3(
    Thermodynamics::SVPFunc1Container &svp_func1,
    std::vector<IndexSet> const &cloud_index_set) {
  Application::Logger app("snap");
  app->Log("Enrolling mgsio3 vapor pressures");

  auto pindex = IndexMap::GetInstance();
  if (!pindex->HasVapor("mgsio3")) return;
  int imgsio3 = pindex->GetVaporId("mgsio3");

  svp_func1[imgsio3][0] = [](AirParcel const &qfrac, int, int) {
    return sat_vapor_p_mgsio3(qfrac.w[IDN]);
  };

  for (int n = 1; n < cloud_index_set[imgsio3].size(); ++n) {
    svp_func1[imgsio3][n] = NullSatVaporPres1;
  }
}

void Thermodynamics::enrollVaporFunctions() {
  enroll_vapor_Function_H2O(svp_func1_, cloud_index_set_);
  enroll_vapor_Function_NH3(svp_func1_, cloud_index_set_);
  enroll_vapor_Function_NH4SH(svp_func1_, cloud_index_set_);
  enroll_vapor_Function_mgsio3(svp_func1_, cloud_index_set_);
}
