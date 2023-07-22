/** @file two_phase_cloud_particles.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Monday May 31, 2021 13:38:00 PDT
 * @bug No known bugs.
 */

// C/C++ headers
#include <iostream>

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "particles.hpp"

TwoPhaseCloudParticles::TwoPhaseCloudParticles(MeshBlock *pmb,
                                               ParameterInput *pin,
                                               std::string name)
    : Particles(pmb, pin, name, 2) {
  ATHENA_LOG("TwoPhaseCloudParticles")
  cnames_.resize(2);
  cnames_[0] = "liquid";
  cnames_[1] = "solid";
  std::cout << "- First category is " << name + " liquid" << std::endl;
  std::cout << "- Second category is " << name + " solid" << std::endl;
}

void TwoPhaseCloudParticles::ExchangeHydro(std::vector<MaterialPoint> &mp,
                                           AthenaArray<Real> &du,
                                           AthenaArray<Real> const &w,
                                           Real dt) {
  Particles::ExchangeHydro(mp, du, w, dt);
  // for (std::vector<MaterialPoint>::iterator it = mp.begin(); it != mp.end();
  // ++it) {
  //   it->v1 += -10.;
  // }
}
