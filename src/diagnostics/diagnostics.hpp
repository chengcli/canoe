#ifndef SRC_DIAGNOSTICS_HPP_
#define SRC_DIAGNOSTICS_HPP_

// C/C++
#include <string>

// Athena++
#include <athena/athena_arrays.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <configure.h>

#include <virtual_groups.hpp>

class ParameterInput;
class Meshlock;

template <typename T, int N>
class PlanarExchanger;

class Diagnostics : public NamedGroup {
 public:
  // data
  std::string type;
  AthenaArray<Real> data;

  // functions
  Diagnostics(MeshBlock *pmb, std::string name);
  virtual ~Diagnostics();

  virtual int GetNumVars() const = 0;
  virtual void Progress(MeshBlock *pmb) {}
  virtual void Finalize(MeshBlock *pmb) {}

 protected:
  int ncells1_, ncells2_, ncells3_;

  //! MPI color of each rank
  std::vector<int> color_;

  //! rank of the bottom block
  std::vector<int> brank_;

  //! scratch geometric arrays
  AthenaArray<Real> x1edge_, x1edge_p1_, x2edge_, x2edge_p1_, x3edge_,
      x3edge_p1_;

  AthenaArray<Real> x1area_, x2area_, x2area_p1_, x3area_, x3area_p1_;

  AthenaArray<Real> vol_;

  std::vector<Real> total_vol_, total_area_;
};

using DiagnosticsPtr = std::shared_ptr<Diagnostics>;
using DiagnosticsContainer = std::vector<DiagnosticsPtr>;

class DiagnosticsFactory {
 public:
  static DiagnosticsContainer CreateFrom(MeshBlock *pmb, ParameterInput *pin);
};

// register all diagnostics
// 1. divergence
class Divergence : public Diagnostics {
 public:
  Divergence(MeshBlock *pmb);
  virtual ~Divergence() {}

  void Finalize(MeshBlock *pmb) override;
  int GetNumVars() const override { return 2; }

 protected:
  AthenaArray<Real> v1f1_, v2f2_, v3f3_;
};

// 2. curl
class Curl : public Diagnostics {
 public:
  Curl(MeshBlock *pmb);
  virtual ~Curl() {}

  void Finalize(MeshBlock *pmb) override;
  int GetNumVars() const override { return 1; }

 protected:
  AthenaArray<Real> v3f2_, v2f3_, v1f3_, v3f1_, v2f1_, v1f2_;
};

// 3. Buoyancy
class Buoyancy : public Diagnostics {
 public:
  Buoyancy(MeshBlock *pmb);
  virtual ~Buoyancy() {}

  void Finalize(MeshBlock *pmb) override;
  int GetNumVars() const override { return 1; }

 protected:
  AthenaArray<Real> pf_;
  Real grav_;
};

// 4. hydro mean
class HydroMean : public Diagnostics {
 public:
  HydroMean(MeshBlock *pmb);
  virtual ~HydroMean() {}

  void Progress(MeshBlock *pmb) override;
  void Finalize(MeshBlock *pmb) override;
  int GetNumVars() const override { return NHYDRO; }

 protected:
  int ncycle_;
};

// 6. anomaly
class Anomaly : public Diagnostics {
 public:
  std::shared_ptr<PlanarExchanger<Real, 2>> pexh;

 public:
  Anomaly(MeshBlock *pmb);
  virtual ~Anomaly() {}

  void Finalize(MeshBlock *pmb) override;
  int GetNumVars() const override { return 4; }

 protected:
  void packData(MeshBlock const *pmb);
  void unpackData(MeshBlock const *pmb);

  //! mean component
  std::vector<Real> mean_;
};

// 8. total radiative flux
/*class RadiativeFlux : public Diagnostics {
 public:
  std::shared_ptr<PlanarExchanger<Real, 2>> pexh;

 public:
  RadiativeFlux(MeshBlock *pmb);
  virtual ~RadiativeFlux() {}

  void Progress(MeshBlock *pmb) override;
  void Finalize(MeshBlock *pmb) override;
  int GetNumVars() const override { return 2; }

 protected:
  void packData(MeshBlock const *pmb);
  void unpackData(MeshBlock const *pmb);

  int ncycle_;
};*/

// 9. hydro flux
class HydroFlux : public Diagnostics {
 public:
  std::shared_ptr<PlanarExchanger<Real, 2>> pexh;

 public:
  HydroFlux(MeshBlock *pmb);
  virtual ~HydroFlux() {}

  void Progress(MeshBlock *pmb) override;
  void Finalize(MeshBlock *pmb) override;
  int GetNumVars() const override { return NHYDRO; }

 protected:
  void packData(MeshBlock const *pmb);
  void unpackData(MeshBlock const *pmb);

  int ncycle_;
};

// 10. w_avg
class V1Moments : public Diagnostics {
 public:
  std::shared_ptr<PlanarExchanger<Real, 2>> pexh;

 public:
  V1Moments(MeshBlock *pmb);
  virtual ~V1Moments() {}

  void Finalize(MeshBlock *pmb) override;
  int GetNumVars() const override { return 3; }

 protected:
  void packData(MeshBlock const *pmb);
  void unpackData(MeshBlock const *pmb);
};

/* 6. eddy flux
class EddyFlux : public Diagnostics {
 public:
  EddyFlux(MeshBlock *pmb);
  virtual ~EddyFlux();
  void Progress(AthenaArray<Real> const &w);
  void Finalize(AthenaArray<Real> const &w);

 protected:
  //! mean and eddy component
  AthenaArray<Real> mean_, eddy_;
};


// 11. total angular momentum
class SphericalAngularMomentum : public Diagnostics {
 public:
  SphericalAngularMomentum(MeshBlock *pmb);
  virtual ~SphericalAngularMomentum() {}
  void Finalize(AthenaArray<Real> const &w);
};

// 12. eddy kinetic energy
class EddyKineticEnergy : public Diagnostics {
 public:
  EddyKineticEnergy(MeshBlock *pmb);
  virtual ~EddyKineticEnergy() {}
  void Finalize(AthenaArray<Real> const &w);
};

// 13. horizontal averaged tendency
class Tendency : public Diagnostics {
 public:
  Tendency(MeshBlock *pmb);
  virtual ~Tendency() {}
  void Finalize(AthenaArray<Real> const &w);

 protected:
  Real last_time_;
  AthenaArray<Real> wh_, up_;
};

class ConvectiveHeatFlux : public Diagnostics {
 public:
  ConvectiveHeatFlux(MeshBlock *pmb);
  virtual ~ConvectiveHeatFlux() {}
  void Progress(AthenaArray<Real> const &w);
  void Finalize(AthenaArray<Real> const &w);
};*/

#endif  // SRC_DIAGNOSTICS_HPP_
