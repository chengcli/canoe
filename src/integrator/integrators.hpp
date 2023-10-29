#ifndef SRC_INTEGRATOR_INTEGRATORS_HPP_
#define SRC_INTEGRATOR_INTEGRATORS_HPP_

// athena
#include <athena/athena.hpp>

struct IntegratorWeight {
  // 2S or 3S* low-storage RK coefficients, Ketchenson (2010)
  Real delta; //!> low-storage coefficients to avoid double F() evaluation per substage
  Real gamma_1, gamma_2, gamma_3; // low-storage coeff for weighted ave of registers
  Real beta; // coeff. from bidiagonal Shu-Osher form Beta matrix, -1 diagonal terms
  Real sbeta, ebeta; // time coeff describing start/end time of each stage
  bool main_stage, orbital_stage; // flag for whether the main calculation is done
};

class MultiStageIntegrator {
 public:
  virtual ~MultiStageIntegrator() {}
  virtual void TimeIntegrate(Real time, Real dt) = 0;
  virtual void WeightedAverage(Real ave_wghts[]) = 0;

  IntegratorWeight const & GetStageWeights(int stage) const {
    return stage_wghts_[stage];
  }

  void SetIntegrator(std::string name);
  int GetNumStages() const { return nstages_; }

 protected:
  IntegratorWeight stage_wghts_[MAX_NSTAGE];
  int nstages_;
};

#endif  // SRC_INTEGRATOR_INTEGRATOR_HPP_
