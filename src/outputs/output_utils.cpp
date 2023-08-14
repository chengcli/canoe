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

  table_ = {// short name, long name, units, grid location
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
            {"radiance", "top-of-atmosphere radiance", "K", "RCC"}};
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
