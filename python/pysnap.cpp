
// snap
#include <snap/thermodynamics/thermodynamics.hpp>

void init_thermo(YAML::Node node) { Thermodynamics::InitFromYAMLInput(node); }

void init_pysnap(py::module &parent) {
  auto m = parent.def_submodule("pysnap", "Python bindings for snap module");

  m.def("init_thermo", &init_thermo);
}
