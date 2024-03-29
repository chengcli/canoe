// pybind
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// snap
#include <snap/thermodynamics/thermodynamics.hpp>

namespace py = pybind11;

void init_snap(py::module &parent) {
  auto m = parent.def_submodule("snap", "Python bindings for snap module");

  m.def("def_thermo", &Thermodynamics::InitFromYAMLInput,
        py::return_value_policy::reference, R"(
      Define thermodynamics for the simulation.

      Parameters
      ----------
      node : YAML::Node
          The thermodynamic configuration node.
      )");

  m.def("def_thermo", &Thermodynamics::InitFromAthenaInput,
        py::return_value_policy::reference);

  py::class_<Thermodynamics>(m, "Thermodynamics")
      .def_static("get", &Thermodynamics::GetInstance)

      .def("get_Rd", &Thermodynamics::GetRd)
      .def("get_GammadRef", &Thermodynamics::GetGammadRef);
}
