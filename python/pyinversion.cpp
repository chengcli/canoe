// pybind
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// harp
#include <harp/radiation.hpp>

// inversion
#include <inversion/profile_inversion.hpp>

namespace py = pybind11;

void init_inversion(py::module &parent) {
  auto m = parent.def_submodule("harp", "Python bindings for inversion module");

  // profile inversion
  py::class_<JunoProfileInversion>(m, "JunoProfileInversion")
      .def(py::init<MeshBlock *, ParameterInput *>())

      .def("cal_fit_target", &JunoProfileInversion::CalculateFitTarget);
}
