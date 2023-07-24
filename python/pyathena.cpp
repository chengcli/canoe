// pybind
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// athena
#include <athena/athena.hpp>
#include <athena/outputs/outputs.hpp>
#include <athena/parameter_input.hpp>

namespace py = pybind11;

PYBIND11_MODULE(pyathena, m) {
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
}
