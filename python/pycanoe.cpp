// pybind
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

// yaml-cpp
#include <yaml-cpp/yaml.h>

// athena
#include <athena/parameter_input.hpp>

// canoe
#include <air_parcel.hpp>
#include <configure.hpp>
#include <constants.hpp>
#include <index_map.hpp>

namespace py = pybind11;

void init_athena(py::module &);
void init_harp(py::module &);
void init_snap(py::module &);
void init_utils(py::module &);

void def_species(std::vector<std::string> vapors,
                 std::vector<std::string> clouds,
                 std::vector<std::string> tracers) {
  IndexMap::InitFromNames(vapors, clouds, tracers);
}

std::string find_resource(const std::string &filename) {
  auto app = Application::GetInstance();
  try {
    auto full_path = app->FindResource(filename);
    return full_path;
  } catch (NotFoundError &e) {
    throw std::runtime_error(e.what());
  }
}

auto cleanup = []() { IndexMap::Destroy(); };

PYBIND11_MODULE(canoe, m) {
  m.attr("__name__") = "canoe";
  m.doc() = "Python bindings for canoe";
  m.add_object("_cleanup", py::capsule(cleanup));

  m.def("def_species", &def_species, py::arg("vapors"),
        py::arg("clouds") = std::vector<std::string>(),
        py::arg("tracers") = std::vector<std::string>());

  m.def("load_configure", &YAML::LoadFile, "");
  m.def("find_resource", &find_resource);

  init_athena(m);
  init_harp(m);
  init_utils(m);
  // init_snap(m);
  //
  //  Constants
  py::module m_constants = m.def_submodule("constants");
  m_constants.attr("Rgas") = Constants::Rgas;
  m_constants.attr("Rgas_cgs") = Constants::Rgas_cgs;
  m_constants.attr("kBoltz") = Constants::kBoltz;
  m_constants.attr("kBoltz_cgs") = Constants::kBoltz_cgs;
  m_constants.attr("Lo") = Constants::Lo;
  m_constants.attr("hPlanck") = Constants::hPlanck;
  m_constants.attr("hPlanck_cgs") = Constants::hPlanck_cgs;
  m_constants.attr("cLight") = Constants::cLight;
  m_constants.attr("cLight_cgs") = Constants::cLight_cgs;
  m_constants.attr("stefanBoltzmann") = Constants::stefanBoltzmann;

  // IndexMap
  py::class_<IndexMap>(m, "index_map")
      .def_static("get", &IndexMap::GetInstance)

      .def("get_vapor_id", &IndexMap::GetVaporId)
      .def("get_cloud_id", &IndexMap::GetCloudId)
      .def("get_tracer_id", &IndexMap::GetTracerId)
      .def("get_species_id", &IndexMap::GetSpeciesId);

  // AirParcel type
  py::enum_<AirParcel::Type>(m, "VariableType")
      .value("MassFrac", AirParcel::Type::MassFrac)
      .value("MassConc", AirParcel::Type::MassConc)
      .value("MoleFrac", AirParcel::Type::MoleFrac)
      .value("MoleConc", AirParcel::Type::MoleConc)
      .export_values();

  // AirParcel
  py::class_<AirParcel>(m, "AirParcel")
      .def(py::init<>())

      .def("hydro",
           [](const AirParcel &var) {
             py::array_t<double> result(NCLOUD, var.c);
             return result;
           })

      .def("cloud",
           [](const AirParcel &var) {
             py::array_t<double> result(NCLOUD, var.c);
             return result;
           })

      .def("tracer",
           [](const AirParcel &var) {
             py::array_t<double> result(NCLOUD, var.x);
             return result;
           })

      .def("to_mass_fraction", &AirParcel::ToMassFraction)
      .def("to_mass_concentration", &AirParcel::ToMassConcentration)
      .def("to_mole_fraction", &AirParcel::ToMoleFraction)
      .def("to_mole_concentration", &AirParcel::ToMoleConcentration);
}
