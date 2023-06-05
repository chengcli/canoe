// C/C++
#include <stdexcept>

// Canoe
#include <configure.hpp>

// athena
#include <athena/mesh/mesh.hpp>
#include <athena/parameter_input.hpp>

// snap
#include <snap/cell_variables.hpp>

#include "mwr_absorbers.hpp"
#include "absorption_functions.hpp"

MwrAbsorberElectron::MwrAbsorberElectron(
  MeshBlock *pmb, ParameterInput *pin, int ion):
    Absorber(pmb, pin, "", "mw_electron")
{
  imols_ = {ion};
  mixrs_ = {1.};
}

Real MwrAbsorberElectron::getAttenuation(Real wave1, Real wave2,
    CellVariables const& var) const
{
  Real P = var.w[IPR]/1.E5; // pa -> bar
  Real T = var.w[IDN];

  Real abs;

  if (model_name_ == "Reference")
    abs = attenuation_freefree_Reference(wave1, P, T);
  else if (model_name_ == "ChengLi")
    abs = attenuation_freefree_Chengli(wave1, P, T);
  else // AppletonHartree
    abs = attenuation_appleton_hartree_nomag(wave1, P, T, var.s[imols_[0]]);

  return 100.*abs;  // 1/cm -> 1/m
}
