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

// communication
#include <communicator/boundary_exchanger.hpp>

// integrator
#include <integrator/integrators.hpp>

class MeshBlock;
class ParameterInput;
class ParticleData;

using ParticleContainer = std::vector<ParticleData>;

class ParticleBase : public NamedGroup,
                     RestartGroup,
                     // ASCIIOutputGroup,
                     // BinaryOutputGroup,
                     MultiStageTimeIntegrator<ParticleContainer>,
                     BoundaryExchanger<ParticleData> {
 public:
  /// public data
  //! particle data container. pc1 is reserved for multi-stage integration
  ParticleContainer pc, pc1;

  /// constructor and destructor
  ParticleBase(MeshBlock *pmb, std::string name);
  virtual ~ParticleBase();

  /// inbound functions
  void SetVelocitiesFromHydro(Hydro const *phydro, Coordinates const *pcoord);

  /// RestartGroup functions
  size_t RestartDataSizeInBytes() const override;
  size_t DumpRestartData(char *pdst) const override;
  size_t LoadRestartData(char *psrt) override;

  /// MultiStageTimeIntegrator functions
  void TimeIntegrate(Real time, Real dt) override;
  void WeightedAverage(ParticleContainer &pc_out,
                       ParticleContainer const &pc_in,
                       Real ave_wghts[]) override;

  /// BoundaryExchanger functions
  void DetachTo(ParticleContainer &buffer) override;
  bool AttachTo(ParticleContainer &container) override;

 protected:
  /// BoundaryExchanger functions
  MeshBlock const *getMeshBlock() const override { return pmy_block_; }

  /// protected data
  //! linked list of particles in cell
  AthenaArray<std::weak_ptr<ParticleData>> pd_in_cell_;

 private:
  //! pointer to parent MeshBlock
  MeshBlock const *pmy_block_;
};

using ParticlePtr = std::shared_ptr<ParticleBase>;
using AllParticles = std::vector<ParticlePtr>;

// helper functions
ParticlePtr find_particle(AllParticles const &pts, std::string name);

// factory methods
namespace ParticlesFactory {
AllParticles create_all_particles(MeshBlock *pmb, ParameterInput *pin);
}  // namespace ParticlesFactory

#endif  // SRC_NBODY_PARTICLE_BASE_HPP_
