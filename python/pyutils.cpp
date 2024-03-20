// pybind
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// yaml-cpp
#include <yaml-cpp/yaml.h>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

namespace py = pybind11;

std::string find_resource(const std::string& filename) {
  auto app = Application::GetInstance();
  try {
    auto full_path = app->FindResource(filename);
    return full_path;
  } catch (NotFoundError& e) {
    throw std::runtime_error(e.what());
  }
}

void init_utils(py::module& parent) {
  auto m = parent.def_submodule("utils", "External utilities");

  m.def("find_resource", &find_resource);

  // yaml-cpp
  py::enum_<YAML::NodeType::value>(m, "NodeType")
      .value("Undefined", YAML::NodeType::Undefined)
      .value("Null", YAML::NodeType::Null)
      .value("Scalar", YAML::NodeType::Scalar)
      .value("Sequence", YAML::NodeType::Sequence)
      .value("Map", YAML::NodeType::Map);

  py::class_<YAML::Node>(m, "YamlNode")
      .def(py::init<const std::string&>())
      .def("__getitem__", [](const YAML::Node node,
                             const std::string& key) { return node[key]; })
      .def("__getitem__",
           [](const YAML::Node node, int n) {
             auto it = node.begin();
             for (int i = 0; i < n; ++i) ++it;
             return *it;
           })
      .def(
          "__iter__",
          [](const YAML::Node& node) {
            return py::make_iterator(node.begin(), node.end());
          },
          py::keep_alive<0, 1>())
      .def("__str__",
           [](const YAML::Node& node) {
             YAML::Emitter out;
             out << node;
             return std::string(out.c_str());
           })
      .def("type", &YAML::Node::Type)
      .def("__len__", &YAML::Node::size);

  py::class_<YAML::detail::iterator_value, YAML::Node>(
      m, "YamlDetailIteratorValue")
      .def(py::init<>())
      .def("first", [](YAML::detail::iterator_value& val) { return val.first; })
      .def("second",
           [](YAML::detail::iterator_value& val) { return val.second; });

  m.def("load_configure", &YAML::LoadFile, "");
}
