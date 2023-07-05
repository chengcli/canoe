// C/C++
#include <cstring>
#include <string>

// athena
#include <athena/parameter_input.hpp>

// application
#include <application/exceptions.hpp>

// canoe
#include <configure.hpp>

void read_thermo_property(Real var[], char const name[], int len, Real v0,
                          ParameterInput* pin) {
  char buf[80], cstr[1024], *p;

  var[0] = v0;
  for (int n = 1; n <= NVAPOR; ++n) {
    snprintf(buf, sizeof(buf), "%s%d", name, n);
    std::string str = pin->GetString("thermodynamics", buf);
    std::snprintf(cstr, sizeof(cstr), "%s", str.c_str());
    p = std::strtok(cstr, " ,");
    int m = 0;
    while ((p != NULL) && (m++ < len)) {
      var[n + (m - 1) * NVAPOR] = std::stod(p);
      p = std::strtok(NULL, " ,");
    }
    if (m != len) {
      throw ValueError("read_thermo_property", name, len, m);
    }
  }
}
