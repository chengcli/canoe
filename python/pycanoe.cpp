// pybind
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// athena
#include <athena/athena.hpp>
#include <athena/parameter_input.hpp>

// canoe
#include <index_map.hpp>

// harp
#include <harp/absorber.hpp>
#include <harp/radiation.hpp>
#include <harp/radiation_band.hpp>

namespace py = pybind11;

PYBIND11_MODULE(pycanoe, m) {
  m.attr("__name__") = "pycanoe";
  m.doc() = "Python bindings for CANOE";

  // IndexMap
  py::class_<IndexMap>(m, "index_map")
      .def_static("get_instance", &IndexMap::GetInstance)
      .def_static("init_from_athena_input", &IndexMap::InitFromAthenaInput)

      .def("get_vapor_id", &IndexMap::GetVaporId)
      .def("get_cloud_id", &IndexMap::GetCloudId)
      .def("get_tracer_id", &IndexMap::GetTracerId)
      .def("get_species_id", &IndexMap::GetSpeciesId);
}
