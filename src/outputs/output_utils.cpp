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
      {"vel", "velocity", "m/s", "CCC"},
      {"vapor", "mass mixing ratio of vapor", "kg/kg", "CCC"},
      {"temp", "temperature", "K", "CCC"},
      {"theta", "potential temperature", "K", "CCC"},
      {"thetav", "virtual potential temperature", "K", "CCC"},
      {"mse", "moist static energy", "J/kg", "CCC"},
      {"rh1", "relative humidity 1", "1", "CCC"},
      {"rh2", "relative humidity 2", "1", "CCC"},
      {"eps", "turbulent dissipation", "w/kg", "CCC"},
      {"tke", "turbulent kinetic energy", "J/kg", "CCC"},
      {"mut", "dynamic turbulent viscosity", "kg/(m.s)", "CCC"},
      {"radiance", "top-of-atmosphere radiance", "K", "RCC"},
      {"curl", "curl", "1/s", "CCC"},
      {"div", "divergence", "1/s", "CCC"},
      {"b", "buoyancy", "m/s^2", "CCC"},
      {"rho_bar", "mean density", "kg/m^3", "CCC"},
      {"q1_bar", "mean vapor1 mixing ratio", "kg/kg", "CCC"},
      {"q2_bar", "mean vapor2 mixing ratio", "kg/kg", "CCC"},
      {"q3_bar", "mean vapor3 mixing ratio", "kg/kg", "CCC"},
      {"vel1_bar", "mean velocity 1", "m/s", "CCC"},
      {"vel2_bar", "mean velocity 2", "m/s", "CCC"},
      {"vel3_bar", "mean velocity 3", "m/s", "CCC"},
      {"T_bar", "mean temperature", "K", "CCC"},
      {"div_h", "horizontal divergence", "1/s", "CCC"},
      {"tempa", "horizontal temperature anomaly", "K", "CCC"},
      {"presa", "horizontal pressure anomaly", "pa", "CCC"},
      {"rflx_up", "total upward radiative flux", "w/m^2", "--F"},
      {"rflx_dn", "total downward radiative flux", "w/m^2", "--F"},
      {"v1rho", "total upward mass flux", "kg/(m^2.s)", "--F"},
      {"v1q1", "total upward vapor1 flux", "1/(m^2.s)", "--C"},
      {"v1q2", "total upward vapor2 flux", "1/(m^2.s)", "--C"},
      {"v1q3", "total upward vapor3 flux", "1/(m^2.s)", "--C"},
      {"v1v1", "total upward velocity 1 flux", "(m/s)/(m^2.s)", "--C"},
      {"v1v2", "total upward velocity 2 flux", "(m/s)/(m^2.s)", "--C"},
      {"v1v3", "total upward velocity 3 flux", "(m/s)/(m^2.s)", "--C"},
      {"v1T", "total upward temperature flux", "K/(m^2.s)", "--C"},
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
