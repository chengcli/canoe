// pybind
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// yaml-cpp
#include <yaml-cpp/yaml.h>

// harp
#include <harp/radiation.hpp>

// inversion
#include <inversion/inversion.hpp>

namespace py = pybind11;

void init_inversion(py::module &parent) {
  auto m = parent.def_submodule("harp", "Python bindings for inversion module");

  // profile inversion
  py::class_<ProfileInversion>(m, "ProfileInversion")
      .def(py::init<YAML::Node const &>())

      .def("update_model", &ProfileInversion::UpdateModel);
}
