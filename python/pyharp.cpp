// pybind
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// athena
#include <athena/athena.hpp>
#include <athena/parameter_input.hpp>

// harp
#include <harp/absorber.hpp>
#include <harp/radiation.hpp>
#include <harp/radiation_band.hpp>

namespace py = pybind11;

PYBIND11_MODULE(pyharp, m) {
  m.attr("__name__") = "pyharp";
  m.doc() = "Python bindings for HARP";

  py::class_<Radiation>(m, "radiation")
      .def_readonly("radiance", &Radiation::radiance)
      .def_readonly("fluxup", &Radiation::flxup)
      .def_readonly("fluxdn", &Radiation::flxdn)

      .def(py::init())

      .def("load_all_radiation_bands", &Radiation::LoadRadiationBands)
      .def("load_radiation_bands", &Radiation::LoadRadiationBands)
      .def("get_num_bands", &Radiation::GetNumBands)
      .def("get_band", &Radiation::GetBand)
      .def("cal_radiative_flux", &Radiation::CalRadiativeFlux)
      .def("cal_radiance", &Radiation::CalRadiance);

  py::class_<RadiationBand>(m, "radband")
      .def_readonly("btau", &RadiationBand::btau)
      .def_readonly("bssa", &RadiationBand::bssa)
      .def_readonly("bpmom", &RadiationBand::bpmom)
      .def_readonly("bflxup", &RadiationBand::bflxup)
      .def_readonly("bflxdn", &RadiationBand::bflxdn)

      .def(py::init())

      .def("get_num_bins", &RadiationBand::GetNumBins)
      .def("has_parameter", &RadiationBand::HasParameter)
      .def("get_parameter", &RadiationBand::GetParameter)
      .def("get_num_absorbers", &RadiationBand::GetNumAbsorbers)
      .def("get_absorber", &RadiationBand::GetAbsorber)
      .def("get_absorber_by_name", &RadiationBand::GetAbsorberByName);
}
