#ifndef SRC_HARP_RT_SOLVERS_HPP_
#define SRC_HARP_RT_SOLVERS_HPP_

// C/C++
#include <string>
#include <Eigen/Dense>  // Needed for Eigen::VectorXd
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <chrono>
#include <utility>

// athena
#include <athena/athena.hpp>

// application
#include <application/application.hpp>
#include <application/exceptions.hpp>

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
  }

  virtual ~RTSolver() {
    Application::Logger app("harp");
    app->Log("Destroy RTSolver " + GetName());
  }

 public:  // member functions
  //! \brief Prepare and seal the solver for the current column
  virtual void Prepare(MeshBlock const *pmb, int k, int j) {}

  //! \brief Allocate memory for radiation solver
  virtual void Resize(int nlyr, int nstr) {
    farea_.DeleteAthenaArray();
    vol_.DeleteAthenaArray();

    farea_.NewAthenaArray(nlyr + 2 * NGHOST);
    vol_.NewAthenaArray(nlyr + 2 * NGHOST);
  }

 public:  // inbound functions
  virtual void CalBandFlux(MeshBlock const *pmb, int k, int j) {
    throw NotImplementedError("CalBandFlux not implemented.");
  }
  virtual void CalBandRadiance(MeshBlock const *pmb, int k, int j) {
    throw NotImplementedError("CalBandRadiance not implemented.");
  }

 protected:
  RadiationBand *pmy_band_;
  AthenaArray<Real> farea_, vol_;
};

class RadiationBand::RTSolverLambert : public RadiationBand::RTSolver {
 public:  // constructor and destructor
  RTSolverLambert(RadiationBand *pmy_band, YAML::Node const &rad)
      : RTSolver(pmy_band, "Lambert") {}
  ~RTSolverLambert() {}

 public:  // inbound functions
  void CalBandRadiance(MeshBlock const *pmb, int k, int j) override;
};

#ifdef RT_DISORT
class RadiationBand::RTSolverDisort : public RadiationBand::RTSolver,
                                      protected DisortWrapper {
 public:  // constructor and destructor
  RTSolverDisort(RadiationBand *pmy_band, YAML::Node const &rad);
  ~RTSolverDisort() {}

 public:  // member functions
  void Prepare(MeshBlock const *pmb, int k, int j) override;
  void Resize(int nlyr, int nstr) override;

 public:  // inbound functions
  void CalBandFlux(MeshBlock const *pmb, int k, int j) override;
  void CalBandRadiance(MeshBlock const *pmb, int k, int j) override;

 protected:
  size_t dir_dim_[2];
  std::vector<Real> dir_axis_;

  void addDisortFlux(Coordinates const *pcoord, int n, int k, int j, int il,
                     int iu);

  void addDisortRadiance(int n, int k, int j);
};


class RadiationBand::RTSolverToon : public RadiationBand::RTSolver,
                                     protected DisortWrapper {
 public:  // Constructor and Destructor
  RTSolverToon(RadiationBand *pmy_band, YAML::Node const &rad);
  ~RTSolverToon() {}

 public:  // Member Functions
  void Prepare(MeshBlock const *pmb, int k, int j) override;
  void Resize(int nlyr, int nstr) override;

 public:  // Inbound Functions
  void CalBandFlux(MeshBlock const *pmb, int k, int j) override;

 protected:
  void addToonFlux(Coordinates const *pcoord, int b, int k, int j, int il, int iu,
                  const Eigen::VectorXd &flux_up, const Eigen::VectorXd &flux_down);
  
  void toonShortwaveSolver(int nlay, double F0_in,
                      const Eigen::VectorXd& mu_in,
                      const Eigen::VectorXd& tau_in,
              	      const Eigen::VectorXd& w_in,
             	      const Eigen::VectorXd& g_in,
             	      double w_surf_in,
                      Eigen::VectorXd& flx_up,
                      Eigen::VectorXd& flx_down);
 
  void toonLongwaveSolver(int nlay,
		      const Eigen::VectorXd& be,
                      const Eigen::VectorXd& tau_in,
                      const Eigen::VectorXd& w_in,
                      const Eigen::VectorXd& g_in,
                      double a_surf_in,
                      Eigen::VectorXd& flx_up,
                      Eigen::VectorXd& flx_down);

  double BB_integrate(double T, double wn1, double wn2);

 private:
      // Private Tridiagonal Solver using the Thomas Algorithm
    inline Eigen::VectorXd tridiagonal_solver(const Eigen::VectorXd& a,  // Sub-diagonal (size n-1)
                                           const Eigen::VectorXd& b,  // Main diagonal (size n)
                                           const Eigen::VectorXd& c,  // Super-diagonal (size n-1)
                                           const Eigen::VectorXd& d); // Right-hand side (size n)
};
  
#endif

#endif  // SRC_HARP_RT_SOLVERS_HPP_
