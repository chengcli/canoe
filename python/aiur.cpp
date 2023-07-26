// pybind
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(aiur, m) {
  py::module::import("pyathena");
  py::module::import("pyharp");
}
