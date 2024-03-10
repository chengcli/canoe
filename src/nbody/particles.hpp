#ifndef SRC_NBODY_PARTICLES_HPP_
#define SRC_NBODY_PARTICLES_HPP_

// C++ headers
#include <memory>
#include <string>
#include <vector>

// Athena++
#include <athena/athena.hpp>

// canoe
#include <virtual_groups.hpp>

// exchanger
#include <exchanger/exchanger.hpp>

// integrator
#include <integrator/integrators.hpp>

class MeshBlock;
class ParameterInput;
class Coordinates;
class ParticleData;
class Hydro;

using ParticleContainer = std::vector<ParticleData>;

class ParticleBase : public NamedGroup,
                     public RestartGroup,
                     // ASCIIOutputGroup,
                     // BinaryOutputGroup,
                     public MeshOutputGroup,
                     public MultiStageIntegrator,
                     public Exchanger<ParticleData, 56> {
 public:
  /// public data
  //! particle data container. pc1 is reserved for multi-stage integration
  ParticleContainer pc, pc1;

  //! mesh data container
  AthenaArray<Real> weight, charge;

 public:  // constructor and destructor
  ParticleBase(MeshBlock *pmb, std::string name);
  virtual ~ParticleBase();

 public:  // member functions
  /// particle-mesh
  void LinkMesh();

 public:  // inbound functions
  void SetVelocitiesFromHydro(Hydro const *phydro, Coordinates const *pcoord);

 public:  // RestartGroup functions
  size_t RestartDataSizeInBytes(Mesh const *pm) const override;
  void DumpRestartData(char *pdst) const override;
  size_t LoadRestartData(char *psrt) override;

 public:  // MeshOutputGroup functions
  bool ShouldMeshOutput(std::string variable_name) const override;
  void LoadMeshOutputData(OutputType *pod, int *num_vars) const override;

 public:  // MultiStageIntegrator functions
  void TimeIntegrate(Real time, Real dt) override;
  void WeightedAverage(Real ave_wghts[]) override;

 public:  // Exchanger functions
  void PackData(MeshBlock const *pmb) override;
  bool UnpackData(MeshBlock const *pmb) override;
  void Transfer(MeshBlock const *pmb, int n = -1) override;
  void ClearBuffer(MeshBlock const *pmb) override;

 protected:
  /// protected data
  //! linked list of particles in cell
  AthenaArray<ParticleData *> pd_in_cell_;

  //! linked flag
  bool linked_flag_;

 private:
  //! pointer to parent MeshBlock
  MeshBlock const *pmy_block_;
};

using ParticlePtr = std::shared_ptr<ParticleBase>;
using AllParticles = std::vector<ParticlePtr>;

// factory methods
class ParticlesFactory {
 public:
  static AllParticles Create(MeshBlock *pmb, ParameterInput *pin);
};

// helper functions
namespace ParticlesHelper {
ParticlePtr find_particle(AllParticles const &pts, std::string name);
}  // namespace ParticlesHelper

namespace AllTasks {

bool launch_particles(MeshBlock *pmb, int stage);
bool integrate_particles(MeshBlock *pmb, int stage);
bool mesh_to_particles(MeshBlock *pmb, int stage);
bool send_particles(MeshBlock *pmb, int stage);
bool recv_particles(MeshBlock *pmb, int stage);
bool particles_to_mesh(MeshBlock *pmb, int stage);
bool attach_particles(MeshBlock *pmb, int stage);

}  // namespace AllTasks

#endif  // SRC_NBODY_PARTICLE_BASE_HPP_
