#ifndef SRC_INVERSION_CONCENTRATION_INVERSION_HPP_
#define SRC_INVERSION_CONCENTRATION_INVERSION_HPP_

// C/C++ header
#include <string>
#include <vector>

// canoe
#include <configure.hpp>

// athena
#include <athena/defs.hpp>

// inversion
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

#endif  // SRC_INVERSION_CONCENTRATION_INVERSION_HPP_
