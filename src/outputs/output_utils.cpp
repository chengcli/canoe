// C/C++
#include <mutex>

// application
#include <application/application.hpp>

// outputs
#include "output_utils.hpp"

static std::mutex table_mutex;

int get_num_variables(std::string grid, AthenaArray<Real> const& data) {
  int nvar;
  if (grid == "--C" || grid == "--F") {
    nvar = data.GetDim2();
  } else if (grid == "---") {
    nvar = data.GetDim1();
  } else {
    nvar = data.GetDim4();
  }

  return nvar;
}

__attribute__((weak)) MetadataTable::MetadataTable() {
  Application::Logger app("outputs");
  app->Log("Initialize MetadataTable");

  table_ = {
      // short name, long name, units, grid location
      {"x1", "height at cell center", "m", "--C"},
      {"x1f", "height at cell boundary", "m", "--F"},
      {"x2", "distance at cell center", "m", "-C-"},
      {"x2f", "distance at cell boundary", "m", "-F-"},
      {"x3", "distance at cell center", "m", "C--"},
      {"x3f", "distance at cell boundary", "m", "F--"},
      {"rho", "density", "kg/m^3", "CCC"},
      {"press", "pressure", "pa", "CCC"},
      {"vel1", "vertical velocity", "m/s", "CCC"},
      {"vel2", "horizontal velocity", "m/s", "CCC"},
      {"vel3", "horizontal velocity", "m/s", "CCC"},
      {"vapor", "mass mixing ratio of vapor", "kg/kg", "CCC"},
      {"temp", "temperature", "K", "CCC"},
      {"theta", "potential temperature", "K", "CCC"},
      {"thetav", "virtual potential temperature", "K", "CCC"},
      {"mse", "moist static energy", "J/kg", "CCC"},
      // relative humidity
      {"rh1", "relative humidity 1", "1", "CCC"},
      {"rh2", "relative humidity 2", "1", "CCC"},
      {"eps", "turbulent dissipation", "w/kg", "CCC"},
      {"tke", "turbulent kinetic energy", "J/kg", "CCC"},
      {"mut", "dynamic turbulent viscosity", "kg/(m.s)", "CCC"},
      // radiation
      {"radiance", "top-of-atmosphere radiance", "K", "RCC"},
      // curl
      {"curl1", "curl in the vertical direction", "1/s", "CCC"},
      {"curl2", "curl in the horizontal direction", "1/s", "CCC"},
      {"curl3", "curl in the horizontal direction", "1/s", "CCC"},
      // divergence
      {"div", "divergence", "1/s", "CCC"},
      {"div_h", "horizontal divergence", "1/s", "CCC"},
      // buoyancy
      {"b", "buoyancy", "m/s^2", "CCC"},
      // hydro mean
      {"rho_bar", "mean density", "kg/m^3", "CCC"},
      {"q1_bar", "mean vapor1 mixing ratio", "kg/kg", "CCC"},
      {"q2_bar", "mean vapor2 mixing ratio", "kg/kg", "CCC"},
      {"q3_bar", "mean vapor3 mixing ratio", "kg/kg", "CCC"},
      {"vel1_bar", "mean velocity 1", "m/s", "CCC"},
      {"vel2_bar", "mean velocity 2", "m/s", "CCC"},
      {"vel3_bar", "mean velocity 3", "m/s", "CCC"},
      {"T_bar", "mean temperature", "K", "CCC"},
      // anomalies
      {"rhoa", "horizontal density anomaly", "kg/m^3", "CCC"},
      {"tempa", "horizontal temperature anomaly", "K", "CCC"},
      {"v1a", "horizontal vertical velocity anomaly", "m/s", "CCC"},
      {"presa", "horizontal pressure anomaly", "pa", "CCC"},
      // radiative flux
      {"rflx_up", "total upward radiative flux", "w/m^2", "--F"},
      {"rflx_dn", "total downward radiative flux", "w/m^2", "--F"},
      // hydro flux
      {"v1rho", "total upward mass flux", "kg/(m^2.s)", "--C"},
      {"v1q1", "total upward vapor1 flux", "1/(m^2.s)", "--C"},
      {"v1q2", "total upward vapor2 flux", "1/(m^2.s)", "--C"},
      {"v1q3", "total upward vapor3 flux", "1/(m^2.s)", "--C"},
      {"v1v1", "total upward velocity 1 flux", "(m/s)/(m^2.s)", "--C"},
      {"v1v2", "total upward velocity 2 flux", "(m/s)/(m^2.s)", "--C"},
      {"v1v3", "total upward velocity 3 flux", "(m/s)/(m^2.s)", "--C"},
      {"v1T", "total upward temperature flux", "K/(m^2.s)", "--C"},
      // v1 moments
      {"w_avg", "mean vertical velocity", "m/s", "--C"},
      {"w2_avg", "mean squard vertical velocity", "(m/s)^2", "--C"},
      {"w3_avg", "mean cubed vertical velocity", "(m/s)^3", "--C"},
      // CMPI6 1D variables description
      {"ta_avg", "domain avg. air temperature profile", "K", "--C"},
      {"ua_avg", "domain avg. zonal wind profile", "m/s", "--C"},
      {"va_avg", "domain avg. meridional wind profile", "m/s", "--C"},
      {"hus_avg", "domain avg. specific humidity profile profile", "kg/kg",
       "--C"},
      {"hur_avg", "domain avg. relative humidity profile", "%", "--C"},
      {"clw_avg", "domain avg. mass fraction of cloud liquid water profile",
       "kg/kg", "--C"},
      {"cli_avg", "domain avg. mass fraction of cloud ice profile", "kg/kg",
       "--C"},
      {"plw_avg",
       "domain avg. mass fraction of precipitating liquid water profile",
       "kg/kg", "--C"},
      {"pli_avg", "domain avg. mass fraction of precipitating ice profile",
       "kg/kg", "--C"},
      {"theta_avg", "domain avg. potential temperature profile", "K", "--C"},
      {"thetae_avg", "domain avg. equivalent potential temperature profile",
       "K", "--C"},
      {"tntrs_avg", "domain avg. shortwave radiative heating rate profile",
       "K/s", "--C"},
      {"tntrl_avg", "domain avg. longwave radiative heating rate profile",
       "K/s", "--C"},
      {"tntrscs_avg",
       "domain avg. shortwave radiative heating rate profile – clear sky",
       "K/s", "--C"},
      {"tntrlcs_avg",
       "domain avg. longwave radiative heating rate profile – clear sky", "K/s",
       "--C"},
      {"cldfrac_avg", "global cloud fraction profile", "%", "--C"},
      // CMPI6 0D variables description
      {"pr_avg", "domain avg. surface precipitation rate", "kg/(m^2 s)", "---"},
      {"hfls_avg", "domain avg. surface upward latent heat flux", "W/m^2",
       "---"},
      {"hfss_avg", "domain avg. surface upward sensible heat flux", "W/m^2",
       "---"},
      {"prw_avg", "domain avg. water vapor path", "kg/m^2", "---"},
      {"clwvi_avg", "domain avg. condensed water path", "kg/m^2", "---"},
      {"clivi_avg", "domain avg. condensed ice water path", "kg/m^2", "---"},
      {"spwr_avg", "domain avg. saturated water vapor path", "kg/m^2", "---"},
      {"rlds_avg", "domain avg. surface downwelling longwave flux", "W/m^2",
       "---"},
      {"rlus_avg", "domain avg. surface upwelling longwave flux", "W/m^2",
       "---"},
      {"rsds_avg", "domain avg. surface downwelling shortwave flux", "W/m^2",
       "---"},
      {"rsus_avg", "domain avg. surface upward shortwave flux", "W/m^2", "---"},
      {"rsdscs_avg",
       "domain avg. surface downwelling shortwave flux – clear sky", "W/m^2",
       "---"},
      {"rsuscs_avg", "domain avg. surface upward shortwave flux – clear sky",
       "W/m^2", "---"},
      {"rldscs_avg",
       "domain avg. surface downwelling longwave flux – clear sky", "W/m^2",
       "---"},
      {"rluscs_avg", "domain avg. surface upwelling longwave flux – clear sky",
       "W/m^2", "---"},
      {"rsdt_avg", "domain avg. TOA incoming shortwave flux", "W/m^2", "---"},
      {"rsut_avg", "domain avg. TOA outgoing shortwave flux", "W/m^2", "---"},
      {"rlut_avg", "domain avg. TOA outgoing longwave flux", "W/m^2", "---"},
      {"rsutcs_avg", "domain avg. TOA outgoing shortwave flux – clear sky",
       "W/m^2", "---"},
      {"rlutcs_avg", "domain avg. TOA outgoing longwave flux – clear sky",
       "W/m^2", "---"},
      // CMPI6 2D variable description
      {"pr", "surface precipitation rate", "kg/(m^2 s)", "CCC"},
      {"pr_conv", "surface convective precipitation rate", "kg/(m^2 s)", "CCC"},
      {"evspsbl", "evaporation flux", "kg/(m^2 s)", "CCC"},
      {"hfls", "surface upward latent heat flux", "W/m^2", "CCC"},
      {"hfss", "surface upward sensible heat flux", "W/m^2", "CCC"},
      {"prw", "water vapor path", "kg/m^2", "CCC"},
      {"clwvi", "condensed water path", "kg/m^2", "CCC"},
      {"clivi", "condensed ice water path", "kg/m^2", "CCC"},
      {"spwr", "saturated water vapor path", "kg/m^2", "CCC"},
      {"rlds", "surface downwelling longwave flux", "W/m^2", "CCC"},
      {"rlus", "surface upwelling longwave flux", "W/m^2", "CCC"},
      {"rsds", "surface downwelling shortwave flux", "W/m^2", "CCC"},
      {"rsus", "surface upward shortwave flux", "W/m^2", "CCC"},
      {"rsdscs", "surface downwelling shortwave flux – clear sky", "W/m^2",
       "CCC"},
      {"rsuscs", "surface upward shortwave flux – clear sky", "W/m^2", "CCC"},
      {"rldscs", "surface downwelling longwave flux – clear sky", "W/m^2",
       "CCC"},
      {"rluscs", "surface upwelling longwave flux – clear sky", "W/m^2", "CCC"},
      {"rsdt", "TOA incoming shortwave flux", "W/m^2", "CCC"},
      {"rsut", "TOA outgoing shortwave flux", "W/m^2", "CCC"},
      {"rlut", "TOA outgoing longwave flux", "W/m^2", "CCC"},
      {"rsutcs", "TOA outgoing shortwave flux – clear sky", "W/m^2", "CCC"},
      {"rlutcs", "TOA outgoing longwave flux – clear sky", "W/m^2", "CCC"},
  };
}

MetadataTable::~MetadataTable() {
  Application::Logger app("outputs");
  app->Log("Destroy MetadataTable");
}

MetadataTable const* MetadataTable::GetInstance() {
  // RAII
  std::unique_lock<std::mutex> lock(table_mutex);

  if (myptr_ == nullptr) {
    myptr_ = new MetadataTable();
  }

  return myptr_;
}

void MetadataTable::Destroy() {
  std::unique_lock<std::mutex> lock(table_mutex);

  if (MetadataTable::myptr_ != nullptr) {
    delete MetadataTable::myptr_;
    MetadataTable::myptr_ = nullptr;
  }
}

std::string MetadataTable::GetGridType(std::string name) const {
  int nouts = table_.size();
  for (int i = 0; i < nouts; ++i) {
    if (table_[i][0] == name) {
      return table_[i][3];
    }
  }

  return "";
}

std::string MetadataTable::GetUnits(std::string name) const {
  int nouts = table_.size();
  for (int i = 0; i < nouts; ++i) {
    if (table_[i][0] == name) {
      return table_[i][2];
    }
  }

  return "";
}

std::string MetadataTable::GetLongName(std::string name) const {
  int nouts = table_.size();
  for (int i = 0; i < nouts; ++i) {
    if (table_[i][0] == name) {
      return table_[i][1];
    }
  }

  return "";
}

MetadataTable* MetadataTable::myptr_ = nullptr;
