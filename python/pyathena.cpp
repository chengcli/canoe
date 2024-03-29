// pybind
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/outputs/outputs.hpp>
#include <athena/parameter_input.hpp>

// canoe
#include <impl.hpp>

// utils
#include <utils/modify_atmoshere.hpp>

namespace py = pybind11;

void init_athena(py::module &parent) {
  auto m = parent.def_submodule("athena", "Python bindings for Athena++");

  m.def("nghost", []() { return NGHOST; });

  py::enum_<IOWrapper::FileMode>(m, "FileMode")
      .value("read", IOWrapper::FileMode::read)
      .value("write", IOWrapper::FileMode::write)
      .export_values();

  py::class_<ParameterInput>(m, "ParameterInput")
      .def(py::init())
      .def("load_from_file",
           [](ParameterInput &pin, const std::string &filename) {
             IOWrapper infile;
             infile.Open(filename.c_str(), IOWrapper::FileMode::read);
             pin.LoadFromFile(infile);
             infile.Close();
           })
      .def("get_integer", &ParameterInput::GetInteger)
      .def("get_real", &ParameterInput::GetReal)
      .def("get_boolean", &ParameterInput::GetBoolean)
      .def("get_string", &ParameterInput::GetString)
      .def("set_string", &ParameterInput::SetString)
      .def("does_parameter_exist", &ParameterInput::DoesParameterExist);

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

  // RegionSize
  py::class_<RegionSize>(m, "RegionSize")
      .def_property(
          "x1min", [](RegionSize const &rs) { return rs.x1min; },
          [](RegionSize &rs, Real x1min) { rs.x1min = x1min; })

      .def_property(
          "x2min", [](RegionSize const &rs) { return rs.x2min; },
          [](RegionSize &rs, Real x2min) { rs.x2min = x2min; })

      .def_property(
          "x3min", [](RegionSize const &rs) { return rs.x3min; },
          [](RegionSize &rs, Real x3min) { rs.x3min = x3min; })

      .def_property(
          "x1max", [](RegionSize const &rs) { return rs.x1max; },
          [](RegionSize &rs, Real x1max) { rs.x1max = x1max; })

      .def_property(
          "x2max", [](RegionSize const &rs) { return rs.x2max; },
          [](RegionSize &rs, Real x2max) { rs.x2max = x2max; })

      .def_property(
          "x3max", [](RegionSize const &rs) { return rs.x3max; },
          [](RegionSize &rs, Real x3max) { rs.x3max = x3max; })

      .def_property(
          "nx1", [](RegionSize const &rs) { return rs.nx1; },
          [](RegionSize &rs, int nx1) { rs.nx1 = nx1; })

      .def_property(
          "nx2", [](RegionSize const &rs) { return rs.nx2; },
          [](RegionSize &rs, int nx2) { rs.nx2 = nx2; })

      .def_property(
          "nx3", [](RegionSize const &rs) { return rs.nx3; },
          [](RegionSize &rs, int nx3) { rs.nx3 = nx3; })

      .def_property(
          "x1rat", [](RegionSize const &rs) { return rs.x1rat; },
          [](RegionSize &rs, Real x1rat) { rs.x1rat = x1rat; })

      .def_property(
          "x2rat", [](RegionSize const &rs) { return rs.x2rat; },
          [](RegionSize &rs, Real x2rat) { rs.x2rat = x2rat; })

      .def_property(
          "x3rat", [](RegionSize const &rs) { return rs.x3rat; },
          [](RegionSize &rs, Real x3rat) { rs.x3rat = x3rat; });

  // Mesh
  py::class_<Mesh>(m, "Mesh")
      .def(py::init<ParameterInput *, int>(), py::arg("pin"),
           py::arg("mesh_only") = false)

      .def("initialize",
           [](Mesh &mesh, ParameterInput *pin) {
             bool restart = false;

             // set up components
             for (int b = 0; b < mesh.nblocal; ++b) {
               MeshBlock *pmb = mesh.my_blocks(b);
               pmb->pimpl = std::make_shared<MeshBlock::Impl>(pmb, pin);
             }
             mesh.Initialize(restart, pin);
           })

      .def("meshblocks",
           [](Mesh &mesh) {
             py::list lst;
             for (size_t i = 0; i < mesh.nbtotal; ++i) {
               lst.append(mesh.my_blocks(i));
             }
             return lst;
           })

      .def(
          "meshblock",
          [](Mesh &mesh, int n) {
            if (n < 0 || n >= mesh.nbtotal) {
              throw py::index_error();
            }
            return mesh.my_blocks(n);
          },
          py::return_value_policy::reference);

  // MeshBlock
  py::class_<MeshBlock>(m, "MeshBlock")
      .def_readonly("block_size", &MeshBlock::block_size)

      .def_readonly("i_st", &MeshBlock::is)
      .def_readonly("i_ed", &MeshBlock::ie)
      .def_readonly("j_st", &MeshBlock::js)
      .def_readonly("j_ed", &MeshBlock::je)
      .def_readonly("k_st", &MeshBlock::ks)
      .def_readonly("k_ed", &MeshBlock::ke)

      //.def_readonly("inversion", [](MeshBlock const& pmb) {
      //  return pmb.pimpl->all_fits;
      //});
      .def("modify_dlnTdlnP",
           [](MeshBlock *mesh_block, Real adlnTdlnP, Real pmin, Real pmax) {
             return modify_atmoshere_adlnTdlnP(mesh_block, adlnTdlnP, pmin,
                                               pmax);
           })

      .def("modify_dlnNH3dlnP",
           [](MeshBlock *mesh_block, Real adlnNH3dlnP, Real pmin, Real pmax) {
             return modify_atmoshere_adlnNH3dlnP(mesh_block, adlnNH3dlnP, pmin,
                                                 pmax);
           });

  // outputs
  py::class_<Outputs>(m, "Outputs")
      .def(py::init<Mesh *, ParameterInput *>(), py::arg("mesh"),
           py::arg("pin"))

      .def("make_outputs", &Outputs::MakeOutputs, py::arg("mesh"),
           py::arg("pin"), py::arg("wtflag") = false);
}
