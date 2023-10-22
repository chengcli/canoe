#ifndef SRC_INTEGRATOR_INTEGRATORS_HPP_
#define SRC_INTEGRATOR_INTEGRATORS_HPP_

template <typename Container>
class MultiStageTimeIntegrator {
 public:
  virtual ~MultiStageTimeIntegrator() {}
  virtual void TimeIntegrate(Real time, Real dt) = 0;
  virtual void WeightedAverage(Container & value_out,
      Container const & value_in,
      Real ave_wghts[]) = 0;
};

#endif  // SRC_INTEGRATOR_INTEGRATOR_HPP_
