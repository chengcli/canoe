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
#include <cppdisort/cppdisort.hpp>
#endif

// harp
#include "radiation_band.hpp"

// communication
#include "column_exchanger.hpp"

class RadiationBand::RTSolver : public NamedGroup {
 public:  // constructor and destructor
  RTSolver(RadiationBand *pmy_band, std::string name)
      : NamedGroup(name), pmy_band_(pmy_band) {
    Application::Logger app("harp");
    app->Log("Initialize RTSolver " + GetName());
  }

  virtual ~RTSolver() {
    Application::Logger app("harp");
    app->Log("Destroy RTSolver " + GetName());
  }

 public:  // inbound functions
  virtual void CalBandFlux(MeshBlock const *pmb, int k, int j, int il, int iu) {
  }

  virtual void CalBandRadiance(MeshBlock const *pmb, int k, int j, int il,
                               int iu) {}

 protected:
  RadiationBand *pmy_band_;
};

class RadiationBand::RTSolverLambert : public RadiationBand::RTSolver {
 public:  // constructor and destructor
  explicit RTSolverLambert(RadiationBand *pmy_band)
      : RTSolver(pmy_band, "Lambert") {}

  ~RTSolverLambert() {}

 public:  // inbound functions
  void CalBandFlux(MeshBlock const *pmb, int k, int j, int il, int iu) override;

  void CalBandRadiance(MeshBlock const *pmb, int k, int j, int il,
                       int iu) override;
};

#ifdef RT_DISORT
class RadiationBand::RTSolverDisort : public RadiationBand::RTSolver,
                                      public DisortWrapper {
 public:  // constructor and destructor
  explicit RTSolverDisort(RadiationBand *pmy_band)
      : RTSolver(pmy_band, "Disort") {}

  ~RTSolverDisort() {}

 public:  // inbound functions
  void CalBandFlux(MeshBlock const *pmb, int k, int j, int il, int iu) override;

  void CalBandRadiance(MeshBlock const *pmb, int k, int j, int il,
                       int iu) override;
};
#endif

#endif  // SRC_HARP_RT_SOLVERS_HPP_
