// pybind
#include <pybind11/numpy.h>
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

void subscribe_species(std::map<std::string, std::vector<std::string>> smap) {
  IndexMap::InitFromSpeciesMap(smap);
}

PYBIND11_MODULE(pyharp, m) {
  m.attr("__name__") = "pyharp";
  m.doc() = "Python bindings for harp module";

  m.def("subscribe_species", &subscribe_species);

  // AthenArray
  py::class_<AthenaArray<Real>>(m, "AthenaArray", py::buffer_protocol())
      .def_buffer([](AthenaArray<Real> &m) -> py::buffer_info {
        size_t stride4 = m.GetDim1() * m.GetDim2() * m.GetDim3() * sizeof(Real);
        size_t stride3 = m.GetDim1() * m.GetDim2() * sizeof(Real);
        size_t stride2 = m.GetDim1() * sizeof(Real);
        size_t stride1 = sizeof(Real);
        if (m.GetDim4() > 1) {
          return py::buffer_info(
              // Pointer to buffer
              m.data(),
              // Size of one scalar
              sizeof(Real),
              // Python struct-style format descriptor
              py::format_descriptor<Real>::format(),
              // Number of dimensions
              4,
              // Buffer dimensions
              {m.GetDim4(), m.GetDim3(), m.GetDim2(), m.GetDim1()},
              // Strides (in bytes) for each index
              {stride4, stride3, stride2, stride1});
        } else if (m.GetDim3() > 1) {
          return py::buffer_info(m.data(), sizeof(Real),
                                 py::format_descriptor<Real>::format(), 3,
                                 {m.GetDim3(), m.GetDim2(), m.GetDim1()},
                                 {stride3, stride2, stride1});
        } else if (m.GetDim2() > 1) {
          return py::buffer_info(
              m.data(), sizeof(Real), py::format_descriptor<Real>::format(), 2,
              {m.GetDim2(), m.GetDim1()}, {stride2, stride1});
        } else {
          return py::buffer_info(m.data(), sizeof(Real),
                                 py::format_descriptor<Real>::format(), 1,
                                 {m.GetDim1()}, {stride1});
        }
      });

  // MeshBlock
  py::class_<MeshBlock>(m, "meshblock");

  // Radiation
  py::class_<Radiation>(m, "radiation")
      .def_readonly("radiance", &Radiation::radiance)
      .def_readonly("fluxup", &Radiation::flxup)
      .def_readonly("fluxdn", &Radiation::flxdn)

      .def(py::init<MeshBlock *, ParameterInput *>())

      .def("get_num_bands", &Radiation::GetNumBands)
      .def("get_band", &Radiation::GetBand)

      .def("cal_flux", &Radiation::CalRadiativeFlux)
      .def("cal_radiance", &Radiation::CalRadiance);

  // RadiationBand
  py::class_<RadiationBand, RadiationBandPtr>(m, "radiation_band")
      .def_readonly("btau", &RadiationBand::btau)
      .def_readonly("bssa", &RadiationBand::bssa)
      .def_readonly("bpmom", &RadiationBand::bpmom)
      .def_readonly("bflxup", &RadiationBand::bflxup)
      .def_readonly("bflxdn", &RadiationBand::bflxdn)

      .def(py::init<std::string, YAML::Node>())
      .def("__str__", &RadiationBand::ToString)

      .def("resize", &RadiationBand::Resize, py::arg("nc1"), py::arg("nc2") = 1,
           py::arg("nc3") = 1, py::arg("nstr") = 4)
      .def("get_num_spec_grids", &RadiationBand::GetNumSpecGrids)
      .def("get_num_absorbers", &RadiationBand::GetNumAbsorbers)
      .def("absorbers", &RadiationBand::Absorbers)
      .def("get_absorber", &RadiationBand::GetAbsorber)
      .def("get_absorber_by_name", &RadiationBand::GetAbsorberByName)
      .def("get_range", &RadiationBand::GetRange)

      .def("cal_flux", &RadiationBand::CalBandFlux, py::arg("pmb") = nullptr,
           py::arg("k") = 0, py::arg("j") = 0, py::arg("is") = 0,
           py::arg("ie") = 0)
      .def("cal_radiance", &RadiationBand::CalBandRadiance,
           py::arg("pmb") = nullptr, py::arg("k") = 0, py::arg("j") = 0)

      .def("get_toa",
           [](RadiationBand &band) {
             py::array_t<Real> ndarray(
                 {band.GetNumSpecGrids(), band.GetNumOutgoingRays()},
                 band._GetToaPtr());
             return ndarray;
           })

      .def("get_tau",
           [](RadiationBand &band) {
             py::array_t<Real> ndarray(
                 {band.GetNumSpecGrids(), band.GetNumLayers()},
                 band._GetTauPtr());
             return ndarray;
           })

      .def("set_spectral_properties",
           [](RadiationBand &band,
              std::map<std::string, std::vector<double>> &atm) {
             auto pindex = IndexMap::GetInstance();

             int nlayer = atm["HGT"].size();

             AirColumn ac(nlayer);
             std::vector<Real> x1v(nlayer), x1f(nlayer + 1);

             for (int i = 0; i < nlayer; ++i) {
               ac[i].SetType(AirParcel::Type::MoleFrac);
               x1v[i] = atm["HGT"][i] * 1e3;  // km -> m
               ac[i].w[IDN] = atm["TEM"][i];
               ac[i].w[IPR] = atm["PRE"][i] * 100.;  // mbar -> Pa
             }

             for (int i = 1; i < nlayer; ++i) {
               x1f[i] = 0.5 * (x1v[i - 1] + x1v[i]);
             }

             x1f[0] = x1v[0] - (x1v[1] - x1v[0]) / 2.;
             x1f[nlayer] =
                 x1v[nlayer - 1] + (x1v[nlayer - 1] - x1v[nlayer - 2]) / 2.;

             for (auto const &pair : atm) {
               if (pair.first == "HGT" || pair.first == "TEM" ||
                   pair.first == "PRE") {
                 continue;
               }

               if (pindex->HasVapor(pair.first))
                 for (int i = 0; i < nlayer; ++i) {
                   ac[i].w[pindex->GetVaporId(pair.first)] =
                       pair.second[i] * 1e-6;
                 }

               if (pindex->HasCloud(pair.first))
                 for (int i = 0; i < nlayer; ++i) {
                   ac[i].c[pindex->GetCloudId(pair.first)] =
                       pair.second[i] * 1e-6;
                 }

               if (pindex->HasChemistry(pair.first))
                 for (int i = 0; i < nlayer; ++i) {
                   ac[i].q[pindex->GetChemistryId(pair.first)] =
                       pair.second[i] * 1e-6;
                 }

               if (pindex->HasTracer(pair.first))
                 for (int i = 0; i < nlayer; ++i) {
                   ac[i].x[pindex->GetTracerId(pair.first)] =
                       pair.second[i] * 1e-6;
                 }
             }

             band.SetSpectralProperties(ac, x1f.data());
           });

  // Absorber
  py::class_<Absorber, AbsorberPtr>(m, "absorber")
      .def("__str__", &Absorber::ToString)
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
           });
}
