//! \file particle_mesh.cpp
//! \brief Implements the functions to transfer particle data to mesh

// C/C++ headers
#include <iostream>

// Athena++ headers
#include <athena/coordinates/coordinates.hpp>
#include <athena/globals.hpp>
#include <athena/mesh/mesh.hpp>

// climath
#include <climath/interpolation.h>  // locate

// n-body
#include "particle_data.hpp"
#include "particles.hpp"

void ParticleBase::LinkMesh() {
  auto pmb = pmy_block_;
  auto pcoord = pmb->pcoord;

  // initialize cell aggregated particle weight u and linked list
  weight.ZeroClear();
  charge.ZeroClear();

  std::fill(pd_in_cell_.data(), pd_in_cell_.data() + pd_in_cell_.GetSize(),
            nullptr);
  int i, j, k;

  auto dims = pcoord->GetDimensions();
  auto x321f = pcoord->GetFaceCoords();

  // loop over particles
  for (auto& qi : pc) {
    k = locate(x321f, qi.x3, dims[0] + 1);
    j = locate(x321f + dims[0] + 1, qi.x2, dims[1] + 1);
    i = locate(x321f + dims[0] + dims[1] + 2, qi.x1, dims[2] + 1);

    auto tid = qi.tid;

    weight(tid, k, j, i) += qi.weight;
    charge(tid, k, j, i) += qi.charge;

    if (pd_in_cell_(tid, k, j, i) == nullptr) {
      pd_in_cell_(tid, k, j, i) = &qi;
      pd_in_cell_(tid, k, j, i)->next = nullptr;
    } else {  // insert sort from low to high weight
      ParticleData* tmp;
      if (pd_in_cell_(tid, k, j, i)->weight < qi.weight) {
        auto qj = pd_in_cell_(tid, k, j, i);
        while (qj->next != nullptr && qj->next->weight < qi.weight)
          qj = qj->next;
        tmp = qj->next;
        qj->next = &qi;
        qj->next->next = tmp;
      } else {
        tmp = pd_in_cell_(tid, k, j, i);
        pd_in_cell_(tid, k, j, i) = &qi;
        pd_in_cell_(tid, k, j, i)->next = tmp;
      }
    }
  }

  linked_flag_ = true;
}
