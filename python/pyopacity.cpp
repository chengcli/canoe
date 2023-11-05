// pybind
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// canoe
#include <air_parcel.hpp>

// opacity
#include <opacity/absorber.hpp>

namespace py = pybind11;

PYBIND11_MODULE(pyopacity, m) {
  m.attr("__name__") = "pyharp";
  m.doc() = "Python bindings for opacity module";
}
