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

  // Absorber
  py::class_<Absorber, AbsorberPtr>(m, "absorber")
      .def("get_attenuation",
           [](Absorber &ab, double wave1, double wave2, py::list lst) {
             AirParcel air;
             std::transform(
                 lst.begin(), lst.end(), air.w,
                 [](const py::handle &elem) { return py::cast<double>(elem); });
             return ab.GetAttenuation(wave1, wave2, air);
           })

      .def("get_single_scattering_albedo",
           [](Absorber &ab, double wave1, double wave2, py::list lst) {
             AirParcel air;
             std::transform(
                 lst.begin(), lst.end(), air.w,
                 [](const py::handle &elem) { return py::cast<double>(elem); });
             return ab.GetSingleScatteringAlbedo(wave1, wave2, air);
           });
}
