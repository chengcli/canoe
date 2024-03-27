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

  // Mesh
  py::class_<Mesh>(m, "Mesh")
      .def(py::init<ParameterInput *, int>(), py::arg("pin"),
           py::arg("mesh_only") = false)

      .def("initialize", [](Mesh &mesh, ParameterInput *pin) {
        bool restart = false;

        // set up components
        for (int b = 0; b < mesh.nblocal; ++b) {
          MeshBlock *pmb = mesh.my_blocks(b);
          pmb->pimpl = std::make_shared<MeshBlock::Impl>(pmb, pin);
        }
        mesh.Initialize(restart, pin);
      });

  // MeshBlock
  py::class_<MeshBlock>(m, "mesh_block");

  // outputs
  py::class_<Outputs>(m, "Outputs")
      .def(py::init<Mesh *, ParameterInput *>(), py::arg("mesh"),
           py::arg("pin"))

      .def("make_outputs", &Outputs::MakeOutputs, py::arg("mesh"),
           py::arg("pin"), py::arg("wtflag") = false);
}
