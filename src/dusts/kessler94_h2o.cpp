// Athena++ header
#include "../mesh/mesh.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "chemistry.hpp"

enum { iH2O = 1, iH2Oc = 2, iH2Op = 3 };

void Chemistry::AssembleReactionMatrix(Chemistry::Vector& r0,
                                       Chemistry::Matrix& r1, Real const q[],
                                       Real time) {
  assert(kc_.size() >= 4);
  Real& k1 = kc_[1];
  Real& k2 = kc_[2];
  Real& k3 = kc_[3];
  r0.setZero();
  r1.setZero();

  Real svp = pmy_block_->pthermo->SatVaporPressure(q[IDN], iH2Oc);
  Real qtol = 1.;  // total gas mols
  for (int n = 1 + NVAPOR; n < ITR; ++n) qtol -= q[n];
  Real dqH2O = q[iH2O] - svp / q[IPR] * qtol;  // q - qs

  if (dqH2O < 0.) {  // evaporation
    r0(iH2O) += -k3 * q[iH2Op] * dqH2O;
    r0(iH2Op) += k3 * q[iH2Op] * dqH2O;
    r1(iH2O, iH2O) += -k3 * q[iH2Op];
    r1(iH2O, iH2Op) += -k3 * dqH2O;
    r1(iH2Op, iH2O) += k3 * q[iH2Op];
    r1(iH2Op, iH2Op) += k3 * dqH2O;
  }

  // autoconversion
  r0(iH2Oc) += -k1 * q[iH2Oc];
  r0(iH2Op) += k1 * q[iH2Oc];
  r1(iH2Oc, iH2Oc) += -k1;
  r1(iH2Op, iH2Oc) += k1;

  // accretion
  r0(iH2Oc) += -k2 * q[iH2Oc] * q[iH2Op];
  r0(iH2Op) += k2 * q[iH2Oc] * q[iH2Op];
  r1(iH2Oc, iH2Oc) += -k2 * q[iH2Op];
  r1(iH2Oc, iH2Op) += -k2 * q[iH2Oc];
  r1(iH2Op, iH2Oc) += k2 * q[iH2Op];
  r1(iH2Op, iH2Op) += k2 * q[iH2Oc];
}
