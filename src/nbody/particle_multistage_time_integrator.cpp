//! \file particle_multistage_time_integrator.cpp
//! \brief Implements functions derived from clas MultistageTimeIntegrator

// C/C++
#include <cstring>

// athena
#include <athena/athena.hpp>

// n-body
#include "particle_data.hpp"
#include "particles.hpp"

void ParticleBase::TimeIntegrate(Real time, Real dt) {
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    for (auto& it : pc) {
      it.x1 += it.v1 * dt;
      it.x2 += it.v2 * dt;
      it.x3 += it.v3 * dt;
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (auto& it : pc) {
      it.x1 += it.v1 * dt;
      it.x2 += it.v2 * dt / it.x1;
      it.x3 += it.v3 * dt / (it.x1 * sin(it.x2));
    }
  }
}

void ParticleBase::WeightedAverage(ParticleContainer& pc_out,
                                   ParticleContainer const& pc_in,
                                   Real ave_wghts[]) {
  size_t psize = pc_out.size();

  for (size_t i = 0; i < psize; ++i) {
    pc_out[i].x1 = ave_wghts[0] * pc_out[i].x1 + ave_wghts[1] * pc_in[i].x1;
    pc_out[i].x2 = ave_wghts[0] * pc_out[i].x2 + ave_wghts[1] * pc_in[i].x2;
    pc_out[i].x3 = ave_wghts[0] * pc_out[i].x3 + ave_wghts[1] * pc_in[i].x3;

    pc_out[i].v1 = ave_wghts[0] * pc_out[i].v1 + ave_wghts[1] * pc_in[i].v1;
    pc_out[i].v2 = ave_wghts[0] * pc_out[i].v2 + ave_wghts[1] * pc_in[i].v2;
    pc_out[i].v3 = ave_wghts[0] * pc_out[i].v3 + ave_wghts[1] * pc_in[i].v3;
  }
}
