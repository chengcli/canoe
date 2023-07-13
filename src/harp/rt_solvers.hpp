#ifndef SRC_HARP_RT_SOLVERS_HPP_
#define SRC_HARP_RT_SOLVERS_HPP_

// C/C++
#include <string>

// athena
#include <athena/athena.hpp>

// application
#include <application/application.hpp>

// canoe
#include <configure.hpp>

// pydisort
#ifdef RT_DISORT
#include <pydisort/cppdisort.hpp>
#endif

// harp
#include "radiation_band.hpp"

class RadiationBand::RTSolver {
 public:
  RTSolver(RadiationBand *pmy_band, std::string name)
      : pmy_band_(pmy_band), name_(name) {
    Application::Logger app("harp");
    app->Log("Initialize RTSolver " + name_);
  }

  std::string GetName() const { return name_; }

  virtual ~RTSolver() {
    Application::Logger app("harp");
    app->Log("Destroy RTSolver " + name_);
  }

  virtual void CalBandFlux(Direction const &rayInput, Real dist_au, int k,
                           int j, int il, int iu) {}

  virtual void CalBandRadiance(Direction const &rayInput, Real dist_au, int k,
                               int j, int il, int iu) {}

 protected:
  RadiationBand *pmy_band_;
  std::string name_;
};

class RadiationBand::RTSolverLambert : public RadiationBand::RTSolver {
 public:
  explicit RTSolverLambert(RadiationBand *pmy_band)
      : RTSolver(pmy_band, "Lambert") {}

  ~RTSolverLambert() {}

  void CalBandFlux(Direction const &rayInput, Real dist_au, int k, int j,
                   int il, int iu) override;

  void CalBandRadiance(Direction const &rayInput, Real dist_au, int k, int j,
                       int il, int iu) override;
};

#ifdef RT_DISORT
class RadiationBand::RTSolverDisort : public RadiationBand::RTSolver,
                                      public DisortWrapper {
 public:
  explicit RTSolverDisort(RadiationBand *pmy_band)
      : RTSolver(pmy_band, "Disort") {}

  ~RTSolverDisort() {}

  void CalBandFlux(Direction const &rayInput, Real dist_au, int k, int j,
                   int il, int iu) override;

  void CalBandRadiance(Direction const &rayInput, Real dist_au, int k, int j,
                       int il, int iu) override;
};
#endif

#endif  // SRC_HARP_RT_SOLVERS_HPP_
