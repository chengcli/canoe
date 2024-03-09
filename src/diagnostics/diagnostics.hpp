#ifndef SRC_DIAGNOSTICS_HPP_
#define SRC_DIAGNOSTICS_HPP_

// C/C++
#include <string>

// Athena++
#include <athena/athena_arrays.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <configure.hpp>
#include <virtual_groups.hpp>

// exchanger
#include <exchanger/exchanger.hpp>

class ParameterInput;

class Diagnostics : public NamedGroup {
 public:
  // data
  std::string type;
  AthenaArray<Real> data;

  // functions
  Diagnostics(MeshBlock *pmb, std::string name);
  virtual ~Diagnostics();

  virtual int GetNumVars() const = 0;
  virtual void Progress(AthenaArray<Real> const &w) {}
  virtual void Finalize(AthenaArray<Real> const &w) {}

 protected:
  MeshBlock *pmy_block_;

  int ncells1_, ncells2_, ncells3_;

  int ncycle_;

  //! mean and eddy component
  AthenaArray<Real> mean_, eddy_;

  //! MPI color of each rank
  std::vector<int> color_;

  //! rank of the bottom block
  std::vector<int> brank_;

  //! scratch geometric arrays
  AthenaArray<Real> x1edge_, x1edge_p1_, x2edge_, x2edge_p1_, x3edge_,
      x3edge_p1_;

  AthenaArray<Real> x1area_, x2area_, x2area_p1_, x3area_, x3area_p1_;

  AthenaArray<Real> vol_, total_vol_;

  void GatherVolumnData(AthenaArray<Real> &total_vol,
                        AthenaArray<Real> &total_data);
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

  void Finalize(AthenaArray<Real> const &w) override;
  int GetNumVars() const override { return 1; }

 protected:
  AthenaArray<Real> v1f1_, v2f2_, v3f3_;
};

// 2. curl
class Curl : public Diagnostics {
 public:
  Curl(MeshBlock *pmb);
  virtual ~Curl() {}

  void Finalize(AthenaArray<Real> const &w) override;
  int GetNumVars() const override { return 1; }

 protected:
  AthenaArray<Real> v3f2_, v2f3_, v1f3_, v3f1_, v2f1_, v1f2_;
};

// 3. Buoyancy
class Buoyancy : public Diagnostics {
 public:
  Buoyancy(MeshBlock *pmb);
  virtual ~Buoyancy() {}

  void Finalize(AthenaArray<Real> const &w) override;
  int GetNumVars() const override { return 1; }

 protected:
  AthenaArray<Real> pf_;
  Real grav_;
};

// 4. hydro mean
class HydroMean : public Diagnostics {
 public:
  HydroMean(MeshBlock *pmb);
  virtual ~HydroMean();

  void Progress(AthenaArray<Real> const &w) override;
  void Finalize(AthenaArray<Real> const &w) override;
  int GetNumVars() const override { return NHYDRO; }
};

// 5. temperature anomaly
class TemperatureAnomaly : public Diagnostics {
 public:
  TemperatureAnomaly(MeshBlock *pmb);
  virtual ~TemperatureAnomaly() {}
  void Finalize(AthenaArray<Real> const &w);
};

// 6. pressure anomaly
class PressureAnomaly : public Diagnostics {
 public:
  PressureAnomaly(MeshBlock *pmb);
  virtual ~PressureAnomaly() {}
  void Finalize(AthenaArray<Real> const &w);

 protected:
  AthenaArray<Real> pf_;
};

/* 6. eddy flux
class EddyFlux : public Diagnostics {
 public:
  EddyFlux(MeshBlock *pmb);
  virtual ~EddyFlux();
  void Progress(AthenaArray<Real> const &w);
  void Finalize(AthenaArray<Real> const &w);
};

// 7. hydro flux
class HydroFlux : public Diagnostics {
 public:
  HydroFlux(MeshBlock *pmb);
  virtual ~HydroFlux() {}
  void Progress(AthenaArray<Real> const &w);
  void Finalize(AthenaArray<Real> const &w);
};

// 8. horizontal divergence
class HorizontalDivergence : public Diagnostics {
 public:
  HorizontalDivergence(MeshBlock *pmb);
  virtual ~HorizontalDivergence();
  void Finalize(AthenaArray<Real> const &w);

 protected:
  AthenaArray<Real> v2f2_, v3f3_;
};


// 10. total radiative flux
class RadiativeFlux : public Diagnostics {
 public:
  RadiativeFlux(MeshBlock *pmb);
  virtual ~RadiativeFlux() {}
  void Progress(AthenaArray<Real> const &w);
  void Finalize(AthenaArray<Real> const &w);
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
