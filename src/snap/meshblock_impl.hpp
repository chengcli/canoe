#ifndef SRC_SNAP_MESHBLOCK_IMPL_HPP_
#define SRC_SNAP_MESHBLOCK_IMPL_HPP_

// C/C++ header
#include <map>
#include <memory>
#include <string>
#include <vector>

// Athena++ header
#include <athena/mesh/mesh.hpp>

class ParameterInput;
class Thermodynamics;
class Decomposition;
class ImplicitSolver;
class FaceReconstruct;
class Forcing;
class Radiation;
class Inversion;

class MeshBlock::IndexMap {
 public:
  IndexMap(MeshBlock *pmb, ParameterInput *pin);
  ~IndexMap() {}

  size_t GetVaporId(std::string name) const {
    return vapor_index_map_.at(name);
  }
  size_t GetTracerId(std::string name) const {
    return tracer_index_map_.at(name);
  }
  size_t GetCloudId(std::string name) const {
    return cloud_index_map_.at(name);
  }
  size_t GetChemistryId(std::string name) const {
    return chemistry_index_map_.at(name);
  }

  size_t GetSpeciesId(std::string category_name) const;

 private:
  MeshBlock *pmy_block_;

  std::map<std::string, size_t> vapor_index_map_;
  std::map<std::string, size_t> tracer_index_map_;
  std::map<std::string, size_t> cloud_index_map_;
  std::map<std::string, size_t> chemistry_index_map_;
  std::map<std::string, size_t> particle_index_map_;
  std::map<std::string, size_t> static_index_map_;
};

//! \class MeshBlock::Impl
//  \brief opaque pointer class implements additional functionality of Athena
//  MeshBlock
class MeshBlock::Impl {
 public:
  Impl(MeshBlock *pmb, ParameterInput *pin);
  ~Impl();

  AthenaArray<Real> du;  // stores tendency

  Thermodynamics *pthermo;
  Decomposition *pdec;
  ImplicitSolver *phevi;
  FaceReconstruct *precon;
  Forcing *pforce;

  Radiation *prad;
  std::vector<Inversion *> fitq;

  Real GetReferencePressure() const { return reference_pressure_; }
  Real GetPressureScaleHeight() const { return pressure_scale_height_; }

 private:
  MeshBlock *pmy_block_;

  Real reference_pressure_;
  Real pressure_scale_height_;
};

#endif  // SRC_SNAP_MESHBLOCK_IMPL_HPP_
