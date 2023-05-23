#ifndef CONCENTRATION_INVERSION_HPP
#define CONCENTRATION_INVERSION_HPP

// C/C++ header
#include <athena/defs.hpp>
#include <configure.hpp>
#include <vector>

#include "inversion.hpp"

class ConcentrationInversion : public Inversion {
 public:
  ConcentrationInversion(MeshBlock *pmb, ParameterInput *pin, std::string name);
  ~ConcentrationInversion();

  void InitializePositions() override;

  void UpdateConcentration(Hydro *phydro, Real *Xp, int k, int jl,
                           int ju) const;

 protected:
  // inversion variable id
  std::vector<int> idx_;

  // prior standard deviation
  Real Xstd_[1 + NVAPOR];
};

#endif
