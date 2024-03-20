// pybind
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// athena
#include <athena/athena.hpp>
#include <athena/mesh/mesh.hpp>
#include <athena/outputs/outputs.hpp>
#include <athena/parameter_input.hpp>

namespace py = pybind11;

void init_athena(py::module &parent) {
  auto m = parent.def_submodule("athena", "Python bindings for Athena++");

  m.def("nghost", []() { return NGHOST; });

  py::enum_<IOWrapper::FileMode>(m, "FileMode")
      .value("read", IOWrapper::FileMode::read)
      .value("write", IOWrapper::FileMode::write)
      .export_values();

  py::class_<IOWrapper>(m, "io_wrapper")
      .def(py::init())
      .def("open", &IOWrapper::Open)
      .def("close", &IOWrapper::Close);

  py::class_<ParameterInput>(m, "parameter_input")
      .def(py::init())
      .def("load_from_file", &ParameterInput::LoadFromFile)
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

  // MeshBlock
  py::class_<MeshBlock>(m, "meshblock");
}
