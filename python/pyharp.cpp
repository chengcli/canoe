// pybind
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// canoe
#include <air_parcel.hpp>
#include <index_map.hpp>

// harp
#include <harp/radiation.hpp>
#include <harp/radiation_band.hpp>

// opacity
#include <opacity/absorber.hpp>

namespace py = pybind11;

void init_index_map(ParameterInput *pin) { IndexMap::InitFromAthenaInput(pin); }

PYBIND11_MODULE(pyharp, m) {
  m.attr("__name__") = "pyharp";
  m.doc() = "Python bindings for harp module";

  m.def("init_index_map", &init_index_map);

  // Radiation
  py::class_<Radiation>(m, "radiation")
      .def_readonly("radiance", &Radiation::radiance)
      .def_readonly("fluxup", &Radiation::flxup)
      .def_readonly("fluxdn", &Radiation::flxdn)

      .def(py::init<MeshBlock *, ParameterInput *>())

      .def("get_num_bands", &Radiation::GetNumBands)
      .def("get_band", &Radiation::GetBand)
      .def("cal_radiative_flux", &Radiation::CalRadiativeFlux)
      .def("cal_radiance", &Radiation::CalRadiance);

  // RadiationBand
  py::class_<RadiationBand, RadiationBandPtr>(m, "radiation_band")
      .def_readonly("btau", &RadiationBand::btau)
      .def_readonly("bssa", &RadiationBand::bssa)
      .def_readonly("bpmom", &RadiationBand::bpmom)
      .def_readonly("bflxup", &RadiationBand::bflxup)
      .def_readonly("bflxdn", &RadiationBand::bflxdn)

      .def(py::init<std::string, YAML::Node>())

      .def("resize", &RadiationBand::Resize)
      .def("get_num_spec_grids", &RadiationBand::GetNumSpecGrids)
      .def("get_num_absorbers", &RadiationBand::GetNumAbsorbers)
      .def("get_absorber", &RadiationBand::GetAbsorber)
      .def("get_absorber_by_name", &RadiationBand::GetAbsorberByName)
      .def("get_range", &RadiationBand::GetRange)

      .def("__str__", [](RadiationBand &band) {
        std::stringstream ss;
        ss << "RadiationBand: " << band.GetName() << std::endl;
        ss << "Absorbers: ";
        for (int n = 0; n < band.GetNumAbsorbers() - 1; ++n) {
          ss << band.GetAbsorber(n)->GetName() << ", ";
        }
        ss << band.GetAbsorber(band.GetNumAbsorbers() - 1)->GetName();
        return ss.str();
      });

  // Absorber
  py::class_<Absorber, AbsorberPtr>(m, "absorber")
      .def("load_opacity_from_file", &Absorber::LoadOpacityFromFile)

      .def("get_name", &Absorber::GetName)

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
           })

      .def("__str__", [](Absorber &ab) {
        std::stringstream ss;
        ss << "Absorber: " << ab.GetName();
        return ss.str();
      });
}
