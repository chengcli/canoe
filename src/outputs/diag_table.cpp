// C/C++
#include <string>
#include <vector>

// outputs
#include "output_utils.hpp"

DiagnosticTable diag_table = {
    // short name, long name, units, grid location
    {"x1", "height at cell center", "m", "--C"},
    {"x1f", "height at cell boundary", "m", "--F"},
    {"x2", "distance at cell center", "m", "-C-"},
    {"x2f", "distance at cell boundary", "m", "-F-"},
    {"x3", "distance at cell center", "m", "C--"},
    {"x3f", "distance at cell boundary", "m", "F--"},
    {"rho", "density", "kg/m^3", "CCC"},
    {"press", "pressure", "pa", "CCC"},
    {"vel", "velocity", "m/s", "CCC"},
    {"vapor", "mass mixing ratio of vapor", "kg/kg", "CCC"},
    {"temp", "temperature", "K", "CCC"},
    {"theta", "potential temperature", "K", "CCC"},
    {"thetav", "virtual potential temperature", "K", "CCC"},
    {"mse", "moist static energy", "J/kg", "CCC"},
    {"rh1", "relative humidity 1", "1", "CCC"},
    {"rh2", "relative humidity 2", "1", "CCC"},
    {"eps", "turbulent dissipation", "w/kg", "CCC"},
    {"tke", "turbulent kinetic energy", "J/kg", "CCC"},
    {"mut", "dynamic turbulent viscosity", "kg/(m.s)", "CCC"},
    {"b1tau", "optical thickness", "1", "CCC"},
    {"b2tau", "optical thickness", "1", "CCC"},
    {"b3tau", "optical thickness", "1", "CCC"},
    {"b4tau", "optical thickness", "1", "CCC"},
    {"b5tau", "optical thickness", "1", "CCC"},
    {"b6tau", "optical thickness", "1", "CCC"},
    {"b1flxup", "upward radiative flux", "w/(m^2)", "CCF"},
    {"b2flxdn", "downward radiative flux", "w/(m^2)", "CCF"},
    {"b3flxup", "upward radiative flux", "w/(m^2)", "CCF"},
    {"b4flxdn", "downward radiative flux", "w/(m^2)", "CCF"},
    {"b5flxup", "upward radiative flux", "w/(m^2)", "CCF"},
    {"b6flxdn", "downward radiative flux", "w/(m^2)", "CCF"},
    {"radiance", "top-of-atmosphere radiance", "K", "RCC"}};
