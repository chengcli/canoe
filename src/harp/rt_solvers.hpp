#ifndef SRC_HARP_RT_SOLVERS_HPP_
#define SRC_HARP_RT_SOLVERS_HPP_

// athena
#include <athena/athena.hpp>

// harp
#include "radiation_band.hpp"

class RadiationBand::RTSolver {
 public:
  RTSolver(RadiationBand *pmy_band);
  ~RTSolver();

  virtual void CalBandFlux(Direction const &rayInput, Real dist_au, int k,
                           int j, int il, int iu) {}

  virtual void CalBandRadiance(Direction const &rayInput, Real dist_au, int k,
                               int j, int il, int iu) {}

 protected:
  RadiationBand *pmy_band_;
};

class RadiationBand::RTSolverLambert : public RadiationBand::RTSolver {
 public:
  RTSolverLambert(RadiationBand *pmy_band);
  ~RTSolverLambert();

  void CalBandFlux(Direction const &rayInput, Real dist_au, int k, int j,
                   int il, int iu) override;

  void CalBandRadiance(Direction const &rayInput, Real dist_au, int k, int j,
                       int il, int iu) override;
};

#endif  // SRC_HARP_RT_SOLVERS_HPP_
