// pybind
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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
void init_utilities(py::module &);

PYBIND11_MODULE(canoe, m) {
  m.attr("__name__") = "canoe";
  m.doc() = "Python bindings for canoe";

  init_athena(m);
  init_harp(m);
  init_utilities(m);
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
      .def_static("from_file", &IndexMap::InitFromAthenaInput)
      .def_static("from_dict", &IndexMap::InitFromSpeciesMap)

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
