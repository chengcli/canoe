// pybind
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// harp
#include <harp/radiation.hpp>

namespace py = pybind11;

PYBIND11_MODULE(xelnaga, m) {
  py::class_<Radiation>(m, "radiation")
      .def_read("radiance", &Radiation::radiance)
      .def_read("fluxup", &Radiation::flxup)
      .def_read("fluxdn", &Radiation::flxdn)

      .def(py::init())

      .def("populate_radiation_bands", &Radiation::PopulateRadiationBands);
}
