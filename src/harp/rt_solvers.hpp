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
#include <virtual_groups.hpp>

// cppdisort
#ifdef RT_DISORT
#include <cppdisort/cppdisort.hpp>
#endif

// harp
#include "radiation_band.hpp"

class RadiationBand::RTSolver : public NamedGroup {
 public:  // constructor and destructor
  RTSolver(RadiationBand *pmy_band, std::string name)
      : NamedGroup(name), pmy_band_(pmy_band) {
    Application::Logger app("harp");
    app->Log("Initialize RTSolver " + GetName());

    farea_.NewAthenaArray(pmy_band_->GetNumLayers() + 2 * NGHOST);
    vol_.NewAthenaArray(pmy_band_->GetNumLayers() + 2 * NGHOST);
  }

  virtual ~RTSolver() {
    Application::Logger app("harp");
    app->Log("Destroy RTSolver " + GetName());
  }

 public:  // member functions
  virtual void Resize(int nlyr, int nstr, int nuphi, int numu) {}

 public:  // inbound functions
  virtual void CalBandFlux(MeshBlock const *pmb, int k, int j, int il, int iu) {
  }

  virtual void CalBandRadiance(MeshBlock const *pmb, int k, int j) {}

 protected:
  RadiationBand *pmy_band_;
  AthenaArray<Real> farea_, vol_;
  int myrank_in_column_;
};

class RadiationBand::RTSolverLambert : public RadiationBand::RTSolver {
 public:  // constructor and destructor
  RTSolverLambert(RadiationBand *pmy_band, YAML::Node const &rad)
      : RTSolver(pmy_band, "Lambert") {}
  ~RTSolverLambert() {}

 public:  // inbound functions
  void CalBandFlux(MeshBlock const *pmb, int k, int j, int il, int iu) override;

  void CalBandRadiance(MeshBlock const *pmb, int k, int j) override;
};

#ifdef RT_DISORT
class RadiationBand::RTSolverDisort : public RadiationBand::RTSolver,
                                      protected DisortWrapper {
 public:  // constructor and destructor
  RTSolverDisort(RadiationBand *pmy_band, YAML::Node const &rad);
  ~RTSolverDisort() {}

 public:  // member functions
  virtual void Resize(int nlyr, int nstr, int nuphi, int numu) override;

 public:  // inbound functions
  void CalBandFlux(MeshBlock const *pmb, int k, int j, int il, int iu) override;

  void CalBandRadiance(MeshBlock const *pmb, int k, int j) override;

 protected:
  void setFlagsFromNode(YAML::Node const &flags);

  Real prepareSpectralBand(MeshBlock const *pmb, int k, int j);

  void addDisortFlux(Coordinates const *pcoord, int n, int k, int j, int il,
                     int iu);

  void addDisortRadiance(Coordinates const *pcoord, int n, int k, int j);
};
#endif

#endif  // SRC_HARP_RT_SOLVERS_HPP_
