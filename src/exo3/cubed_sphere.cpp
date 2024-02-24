// C/C++
#include <cmath>
#include <iostream>
#include <sstream>

// athena
#include <athena/coordinates/coordinates.hpp>
#include <athena/globals.hpp>
#include <athena/hydro/hydro.hpp>
#include <athena/mesh/mesh.hpp>

// canoe
#include <configure.hpp>

// exo3
#include "cubed_sphere.hpp"
#include "cubed_sphere_utility.hpp"

// check
#include <checks.hpp>

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

namespace cs = CubedSphereUtility;

CubedSphere::CubedSphere(MeshBlock *pmb) : pmy_block_(pmb) {
  int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;

  L3DValues[0].NewAthenaArray(NWAVE, nc3, nc2, nc1);
  L3DValues[1].NewAthenaArray(NWAVE, nc3, nc2, nc1);
  L3DValues[2].NewAthenaArray(NWAVE, nc3, nc2, nc1);

  R3DValues[0].NewAthenaArray(NWAVE, nc3, nc2, nc1);
  R3DValues[1].NewAthenaArray(NWAVE, nc3, nc2, nc1);
  R3DValues[2].NewAthenaArray(NWAVE, nc3, nc2, nc1);
}

Real CubedSphere::GenerateMeshX2(Real x, LogicalLocation const &loc) {
  Real x_l, x_u;
  int lx2_lv2 = loc.lx2 >> (loc.level - 2);

  switch (lx2_lv2) {
    case 0:
      x_l = -0.5;
      x_u = 0.0;
      break;
    case 1:
      x_l = 0.0;
      x_u = 0.5;
  }

  return (0.5 * (x - x_l) / (x_u - x_l) - 0.25) * PI;  // Add Pi back later!!
}

Real CubedSphere::GenerateMeshX3(Real x, LogicalLocation const &loc) {
  Real x_l, x_u;
  int lx3_lv2 = loc.lx3 >> (loc.level - 2);

  switch (lx3_lv2) {
    case 0:
      x_l = -0.5;
      x_u = -1.0 / 6.0;
      break;
    case 1:
      x_l = -1.0 / 6.0;
      x_u = 1.0 / 6.0;
      break;
    case 2:
      x_l = 1.0 / 6.0;
      x_u = 0.5;
  }

  return (0.5 * (x - x_l) / (x_u - x_l) - 0.25) * PI;  // Add Pi back later!!
}

// Obtain Lat and Lon (radians) from x2 and x3
// k is not used for now
// Find the block number
void CubedSphere::GetLatLon(Real *lat, Real *lon, int k, int j, int i) const {
  auto pcoord = pmy_block_->pcoord;
  auto &loc = pmy_block_->loc;

  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int blockID = lv2_lx2 + lv2_lx3 * 2 + 1;
  // Calculate the needed parameters
  Real dX = tan(pcoord->x2v(j));
  Real dY = tan(pcoord->x3v(k));
  cs::RLLFromXYP(dY, -dX, blockID - 1, *lon, *lat);
  *lon = 2.0 * PI - *lon;
}

// Obtain Lat and Lon (radians) from x2 and x3
// k is not used for now
// Find the block number
void CubedSphere::GetLatLonFace2(Real *lat, Real *lon, int k, int j,
                                 int i) const {
  auto pcoord = pmy_block_->pcoord;
  auto &loc = pmy_block_->loc;

  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int blockID = lv2_lx2 + lv2_lx3 * 2 + 1;

  // Calculate the needed parameters
  Real dX = tan(pcoord->x2f(j));
  Real dY = tan(pcoord->x3v(k));
  cs::RLLFromXYP(dY, -dX, blockID - 1, *lon, *lat);
  *lon = 2.0 * PI - *lon;
}

// Obtain Lat and Lon (radians) from x2 and x3
// k is not used for now
// Find the block number
void CubedSphere::GetLatLonFace3(Real *lat, Real *lon, int k, int j,
                                 int i) const {
  auto pcoord = pmy_block_->pcoord;
  auto &loc = pcoord->pmy_block->loc;

  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int blockID = lv2_lx2 + lv2_lx3 * 2 + 1;

  // Calculate the needed parameters
  Real dX = tan(pcoord->x2v(j));
  Real dY = tan(pcoord->x3f(k));
  cs::RLLFromXYP(dY, -dX, blockID - 1, *lon, *lat);
  *lon = 2.0 * PI - *lon;
}

// Obtain U and V (Lat-Lon) from V2 and V3 (Gnomonic Equiangle)
// U is V_lam, V is V_phi.
void CubedSphere::GetUV(Real *U, Real *V, Real V2, Real V3, int k, int j,
                        int i) const {
  auto pcoord = pmy_block_->pcoord;
  auto &loc = pmy_block_->loc;

  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int blockID = lv2_lx2 + lv2_lx3 * 2 + 1;

  Real X = tan(pcoord->x2v(j));
  Real Y = tan(pcoord->x3v(k));
  cs::VecTransRLLFromABP(X, Y, blockID, V2, V3, U, V);
  *U = -*U;
}

// Convert U and V (Lat-Lon) to V2 and V3 (Gnomonic Equiangle)
// U is V_lam, V is V_phi.
void CubedSphere::GetVyVz(Real *V2, Real *V3, Real U, Real V, int k, int j,
                          int i) const {
  auto pcoord = pmy_block_->pcoord;
  auto &loc = pmy_block_->loc;

  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int blockID = lv2_lx2 + lv2_lx3 * 2 + 1;

  // Calculate the needed parameters
  Real X = tan(pcoord->x2v(j));
  Real Y = tan(pcoord->x3v(k));
  cs::VecTransABPFromRLL(X, Y, blockID, -U, V, V2, V3);
}

void CubedSphere::CalculateCoriolisForce2(int i2, int i3, Real v2, Real v3,
                                          Real Omega, Real den, Real *cF2,
                                          Real *cF3) const {
  auto pcoord = pmy_block_->pcoord;
  auto &loc = pmy_block_->loc;

  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int blockID = lv2_lx2 + lv2_lx3 * 2 + 1;

  Real x = tan(pcoord->x2v(i2));
  Real y = tan(pcoord->x3v(i3));
  Real C = sqrt(x * x + 1);
  Real D = sqrt(y * y + 1);
  Real delta = sqrt(x * x + y * y + 1);
  Real A = 2 * Omega * C * D / delta / delta * den;

  switch (blockID) {
    case 1:
    case 5:
      *cF2 = -A * v3;
      *cF3 = A * v2;
      break;
    case 2:
    case 3:
    case 4:
    case 6:
      *cF2 = A * x * v3;
      *cF3 = -x * A * v2;
      break;
  }
}

void CubedSphere::CalculateCoriolisForce3(int i2, int i3, Real v1, Real v2,
                                          Real v3, Real Omega, Real den,
                                          Real *cF1, Real *cF2,
                                          Real *cF3) const {
  auto pcoord = pmy_block_->pcoord;
  auto &loc = pmy_block_->loc;

  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int blockID = lv2_lx2 + lv2_lx3 * 2 + 1;

  Real x = tan(pcoord->x2v(i2));
  Real y = tan(pcoord->x3v(i3));
  Real C = sqrt(x * x + 1);
  Real D = sqrt(y * y + 1);
  Real delta = sqrt(x * x + y * y + 1);
  Real A = 2 * Omega * C * D / delta / delta * den;

  switch (blockID) {
    case 1:
    case 5:
      *cF1 = A * (y * C * v2 - x * D * v3);
      *cF2 = A * (-y * C * v1 - v3);
      *cF3 = A * (x * D * v1 + v2);
      break;
    case 2:
    case 3:
    case 4:
    case 6:
      *cF1 = -A * D * v3;
      *cF2 = A * x * v3;
      *cF3 = A * (D * v1 - x * v2);
      break;
  }
}

void CubedSphere::CovariantVectorToContravariant(int i2, int i3, Real v2,
                                                 Real v3, Real *v2c,
                                                 Real *v3c) const {
  auto pcoord = pmy_block_->pcoord;

  Real x = tan(pcoord->x2v(i2));
  Real y = tan(pcoord->x3v(i3));
  Real C = sqrt(x * x + 1);
  Real D = sqrt(y * y + 1);
  Real delta = sqrt(x * x + y * y + 1);
  Real cth = -x * y / (C * D);
  Real sth = delta / (C * D);

  *v2c = v2 / sth / sth - v3 * cth / sth / sth;
  *v3c = v3 / sth / sth - v2 * cth / sth / sth;
}

void CubedSphere::ContravariantVectorToCovariant(int i2, int i3, Real v2,
                                                 Real v3, Real *v2c,
                                                 Real *v3c) const {
  auto pcoord = pmy_block_->pcoord;

  Real x = tan(pcoord->x2v(i2));
  Real y = tan(pcoord->x3v(i3));
  Real C = sqrt(x * x + 1);
  Real D = sqrt(y * y + 1);
  Real delta = sqrt(x * x + y * y + 1);
  Real cth = -x * y / (C * D);
  Real sth = delta / (C * D);

  *v2c = v2 + v3 * cth;
  *v3c = v3 + v2 * cth;
}

void CubedSphere::SaveLR3DValues(AthenaArray<Real> &L_in,
                                 AthenaArray<Real> &R_in, int direction, int k,
                                 int j, int il, int iu) {
  for (int n = 0; n < NWAVE; n++) {
    for (int i = il; i <= iu; i++) {
      L3DValues[direction](n, k, j, i) = L_in(n, i);
      R3DValues[direction](n, k, j, i) = R_in(n, i);
    }
  }
}

void CubedSphere::LoadLR3DValues(AthenaArray<Real> &L_in,
                                 AthenaArray<Real> &R_in, int direction, int k,
                                 int j, int il, int iu) {
  for (int n = 0; n < NWAVE; n++) {
    for (int i = il; i <= iu; i++) {
      L_in(n, i) = L3DValues[direction](n, k, j, i);
      R_in(n, i) = R3DValues[direction](n, k, j, i);
    }
  }
}

void CubedSphere::SynchronizeFluxesSend() {
#ifdef MPI_PARALLEL
  MeshBlock *pmb = pmy_block_;
  for (int i = 0; i < 4; ++i) send_flag_[i] = 0;

  for (int n = 0; n < pmb->pbval->nneighbor; n++) {
    NeighborBlock &nb = pmb->pbval->neighbor[n];
    if (nb.ni.ox1 == 0 &&
        nb.ni.ox2 * nb.ni.ox3 == 0) {  // On x2 and x3 face boundaries only
      sendNeighborBlocks(nb.ni.ox2, nb.ni.ox3, nb.snb.rank, nb.snb.gid);
    }
  }
#endif
}

void CubedSphere::SynchronizeFluxesRecv() {
#ifdef MPI_PARALLEL
  MeshBlock *pmb = pmy_block_;
  for (int i = 0; i < 4; ++i) recv_flag_[i] = 0;

  for (int n = 0; n < pmb->pbval->nneighbor; n++) {
    NeighborBlock &nb = pmb->pbval->neighbor[n];
    if (nb.ni.ox1 == 0 &&
        nb.ni.ox2 * nb.ni.ox3 == 0) {  // On x2 and x3 face boundaries only
      recvNeighborBlocks(nb.ni.ox2, nb.ni.ox3, nb.snb.rank, nb.snb.gid);
    }
  }
#endif
}

void CubedSphere::SynchronizeFluxesWait() {
#ifdef MPI_PARALLEL
  MPI_Status status;

  for (int i = 0; i < 4; ++i) {
    if (send_flag_[i]) MPI_Wait(&send_request_[i], &status);
  }
#endif
}

void CubedSphere::sendNeighborBlocks(int ox2, int ox3, int tg_rank,
                                     int tg_gid) {
  MeshBlock *pmb = pmy_block_;
  auto &loc = pmb->loc;

  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int blockID = lv2_lx2 + lv2_lx3 * 2 + 1;
  // Calculate local ID
  int local_lx2 = loc.lx2 - (lv2_lx2 << (loc.level - 2));
  int local_lx3 = loc.lx3 - (lv2_lx3 << (loc.level - 2));
  int bound_lim = (1 << (loc.level - 2)) - 1;
  int tox2, tox3;
  int DirTag, DirNum, ownTag;  // Tag for target, and the numbering of axis
  int ox2_bkp = ox2;
  int ox3_bkp = ox3;
  bool Left;       // Marking whether need to reverse direction in packing or
                   // unpacking; Left marking left or right.
  int target_dir;  // Up=0, Down=1, Left=2, Right=3

  // Target block ID
  const int tgbid[6][4] = {// To access: tgbid[source_id][target_dir]
                           {6, 2, 3, 4}, {1, 5, 3, 4}, {1, 5, 6, 2},
                           {1, 5, 2, 6}, {2, 6, 3, 4}, {1, 5, 4, 3}};

  // Direction inverse table
  const int dinv[6][4] = {// To access: dinv[source_id][target_dir]
                          {1, 0, 0, 1}, {0, 0, 0, 0}, {0, 1, 0, 0},
                          {1, 0, 0, 0}, {0, 1, 1, 0}, {1, 1, 0, 0}};

  // Hard code the boundary cases
  if (local_lx2 == bound_lim && ox2 == 1) {  // Down
    target_dir = 1;
    TransformOX(&ox2_bkp, &ox3_bkp, &tox2, &tox3);
    Left = true;
  } else if (local_lx2 == 0 && ox2 == -1) {  // Up
    target_dir = 0;
    TransformOX(&ox2_bkp, &ox3_bkp, &tox2, &tox3);
    Left = false;
  } else if (local_lx3 == bound_lim && ox3 == 1) {  // Right
    target_dir = 3;
    TransformOX(&ox2_bkp, &ox3_bkp, &tox2, &tox3);
    Left = true;
  } else if (local_lx3 == 0 && ox3 == -1) {  // Left
    target_dir = 2;
    TransformOX(&ox2_bkp, &ox3_bkp, &tox2, &tox3);
    Left = false;
  } else {  // No need to communicate fluxes, return
    return;
  }
  // Pack the data
  int kb1, kb2, jb1, jb2, ib1, ib2;
  if (ox2 == 1) {
    DirNum = X2DIR;
    jb1 = pmb->je + 1;
    jb2 = pmb->je + 1;
    ib1 = pmb->is - 1;
    ib2 = pmb->ie + 1;
    kb1 = pmb->ks;
    kb2 = pmb->ke;
  }
  if (ox2 == -1) {
    DirNum = X2DIR;
    jb1 = pmb->js;
    jb2 = pmb->js;
    ib1 = pmb->is - 1;
    ib2 = pmb->ie + 1;
    kb1 = pmb->ks;
    kb2 = pmb->ke;
  }
  if (ox3 == 1) {
    DirNum = X3DIR;
    jb1 = pmb->js;
    jb2 = pmb->je;
    ib1 = pmb->is - 1;
    ib2 = pmb->ie + 1;
    kb1 = pmb->ke + 1;
    kb2 = pmb->ke + 1;
  }
  if (ox3 == -1) {
    DirNum = X3DIR;
    jb1 = pmb->js;
    jb2 = pmb->je;
    ib1 = pmb->is - 1;
    ib2 = pmb->ie + 1;
    kb1 = pmb->ks;
    kb2 = pmb->ks;
  }
  int dsize = ((kb2 - kb1 + 1) * (jb2 - jb1 + 1) * (ib2 - ib1 + 1) * NWAVE);
  LRDataBuffer[target_dir].resize(dsize);
  Real *data = LRDataBuffer[target_dir].data();
  int offset = 0;
  // Added Feb 25: Calculate the basis vectors of source & target panels
  if (dinv[blockID - 1][target_dir] == 1) {
    for (int n = 0; n < NWAVE; n++)
      for (int k = kb2; k >= kb1; k--)
        for (int j = jb2; j >= jb1; j--)
          for (int i = ib1; i <= ib2; i++)
            if ((n == IVY) || (n == IVZ)) {  // Projection needed
              Real x2s, x3s, x2d, x3d;       // Coordinates of source and dest
              Real pos_along;
              if (ox2 == 1) {
                x2s = PI / 4.0;
                x3s = pmb->pcoord->x3v(
                    k);  // Need to check, is k starting in the right position?
                pos_along = x3s;
              } else if (ox2 == -1) {
                x2s = -PI / 4.0;
                x3s = pmb->pcoord->x3v(k);
                pos_along = x3s;
              } else if (ox3 == 1) {
                x2s = pmb->pcoord->x2v(j);
                x3s = PI / 4.0;
                pos_along = x2s;
              } else if (ox3 == -1) {
                x2s = pmb->pcoord->x2v(j);
                x3s = -PI / 4.0;
                pos_along = x2s;
              }
              if (tox2 == 1) {
                x2d = PI / 4.0;
                x3d = -pos_along;
              } else if (tox2 == -1) {
                x2d = -PI / 4.0;
                x3d = -pos_along;
              } else if (tox3 == 1) {
                x2d = -pos_along;
                x3d = PI / 4.0;
              } else if (tox3 == -1) {
                x2d = -pos_along;
                x3d = -PI / 4.0;
              }
              x2d = tan(x2d);
              x3d = tan(x3d);
              x2s = tan(x2s);
              x3s = tan(x3s);
              Real U, V, vx, vy;
              if (Left) {
                cs::VecTransRLLFromABP(x2s, x3s, blockID,
                                       L3DValues[DirNum](IVY, k, j, i),
                                       L3DValues[DirNum](IVZ, k, j, i), &U, &V);
                cs::VecTransABPFromRLL(x2d, x3d, tgbid[blockID - 1][target_dir],
                                       U, V, &vx, &vy);
              } else {
                cs::VecTransRLLFromABP(x2s, x3s, blockID,
                                       R3DValues[DirNum](IVY, k, j, i),
                                       R3DValues[DirNum](IVZ, k, j, i), &U, &V);
                cs::VecTransABPFromRLL(x2d, x3d, tgbid[blockID - 1][target_dir],
                                       U, V, &vx, &vy);
              }
              if (n == IVY)
                data[offset++] = vx;
              else
                data[offset++] = vy;
            } else {  // No projection needed
              if (Left)
                data[offset++] = L3DValues[DirNum](n, k, j, i);
              else
                data[offset++] = R3DValues[DirNum](n, k, j, i);
            }
  } else {
    for (int n = 0; n < NWAVE; n++)
      for (int k = kb1; k <= kb2; k++)
        for (int j = jb1; j <= jb2; j++)
          for (int i = ib1; i <= ib2; i++)
            if ((n == IVY) || (n == IVZ)) {  // Projection needed
              Real x2s, x3s, x2d, x3d;       // Coordinates of source and dest
              Real pos_along;
              if (ox2 == 1) {
                x2s = PI / 4.0;
                x3s = pmb->pcoord->x3v(
                    k);  // Need to check, is k starting in the right position?
                pos_along = x3s;
              } else if (ox2 == -1) {
                x2s = -PI / 4.0;
                x3s = pmb->pcoord->x3v(k);
                pos_along = x3s;
              } else if (ox3 == 1) {
                x2s = pmb->pcoord->x2v(j);
                x3s = PI / 4.0;
                pos_along = x2s;
              } else if (ox3 == -1) {
                x2s = pmb->pcoord->x2v(j);
                x3s = -PI / 4.0;
                pos_along = x2s;
              }
              if (tox2 == 1) {
                x2d = PI / 4.0;
                x3d = pos_along;
              } else if (tox2 == -1) {
                x2d = -PI / 4.0;
                x3d = pos_along;
              } else if (tox3 == 1) {
                x2d = pos_along;
                x3d = PI / 4.0;
              } else if (tox3 == -1) {
                x2d = pos_along;
                x3d = -PI / 4.0;
              }
              // Angular positions to tan values
              x2d = tan(x2d);
              x2s = tan(x2s);
              x3d = tan(x3d);
              x3s = tan(x3s);
              Real U, V, vx, vy;
              if (Left) {
                cs::VecTransRLLFromABP(x2s, x3s, blockID,
                                       L3DValues[DirNum](IVY, k, j, i),
                                       L3DValues[DirNum](IVZ, k, j, i), &U, &V);
                cs::VecTransABPFromRLL(x2d, x3d, tgbid[blockID - 1][target_dir],
                                       U, V, &vx, &vy);
              } else {
                cs::VecTransRLLFromABP(x2s, x3s, blockID,
                                       R3DValues[DirNum](IVY, k, j, i),
                                       R3DValues[DirNum](IVZ, k, j, i), &U, &V);
                cs::VecTransABPFromRLL(x2d, x3d, tgbid[blockID - 1][target_dir],
                                       U, V, &vx, &vy);
              }
              if (n == IVY)
                data[offset++] = vx;
              else
                data[offset++] = vy;
            } else {  // No projection needed
              if (Left)
                data[offset++] = L3DValues[DirNum](n, k, j, i);
              else
                data[offset++] = R3DValues[DirNum](n, k, j, i);
            }
  }

  // Calculate the tag of destination
  if (tox2 == -1)
    DirTag = 0 + 4 * pmb->gid + 24 * (1 << (loc.level - 2)) * tg_gid;
  if (tox2 == 1)
    DirTag = 1 + 4 * pmb->gid + 24 * (1 << (loc.level - 2)) * tg_gid;
  if (tox3 == -1)
    DirTag = 2 + 4 * pmb->gid + 24 * (1 << (loc.level - 2)) * tg_gid;
  if (tox3 == 1)
    DirTag = 3 + 4 * pmb->gid + 24 * (1 << (loc.level - 2)) * tg_gid;
  // Send by MPI: we don't care whether it is in the same process for now
  if (ox2 == -1) ownTag = 0;
  if (ox2 == 1) ownTag = 1;
  if (ox3 == -1) ownTag = 2;
  if (ox3 == 1) ownTag = 3;

#ifdef MPI_PARALLEL
  MPI_Isend(data, dsize, MPI_DOUBLE, tg_rank, DirTag, MPI_COMM_WORLD,
            &send_request_[ownTag]);
  send_flag_[ownTag] = 1;
#endif
}

void CubedSphere::recvNeighborBlocks(int ox2, int ox3, int tg_rank,
                                     int tg_gid) {
  MeshBlock *pmb = pmy_block_;
  auto &loc = pmb->loc;

  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int blockID = FindBlockID(loc);
  // Calculate local ID
  int local_lx2 = loc.lx2 - (lv2_lx2 << (loc.level - 2));
  int local_lx3 = loc.lx3 - (lv2_lx3 << (loc.level - 2));
  int bound_lim = (1 << (loc.level - 2)) - 1;
  int tox2, tox3;
  int DirTag, DirNum, ownTag;  // Tag for receiving, and the numbering of axis
                               // to place the values
  int ox2_bkp = ox2;
  int ox3_bkp = ox3;
  bool invDir, Left;  // Marking whether need to reverse direction in packing or
                      // unpacking; Left marking left or right.
  // Hard code the boundary cases
  if (local_lx2 == bound_lim && ox2 == 1) {
    TransformOX(&ox2_bkp, &ox3_bkp, &tox2, &tox3);
    Left = false;
    // Determine whether need to inverse direction
  } else if (local_lx2 == 0 && ox2 == -1) {
    TransformOX(&ox2_bkp, &ox3_bkp, &tox2, &tox3);
    Left = true;
  } else if (local_lx3 == bound_lim && ox3 == 1) {
    TransformOX(&ox2_bkp, &ox3_bkp, &tox2, &tox3);
    Left = false;
  } else if (local_lx3 == 0 && ox3 == -1) {
    TransformOX(&ox2_bkp, &ox3_bkp, &tox2, &tox3);
    Left = true;
  } else {  // No need to communicate fluxes, return
    return;
  }
  // Pack the data
  int kb1, kb2, jb1, jb2, ib1, ib2;
  if (ox2 == 1) {
    DirNum = X2DIR;
    jb1 = pmb->je + 1;
    jb2 = pmb->je + 1;
    ib1 = pmb->is - 1;
    ib2 = pmb->ie + 1;
    kb1 = pmb->ks;
    kb2 = pmb->ke;
  }
  if (ox2 == -1) {
    DirNum = X2DIR;
    jb1 = pmb->js;
    jb2 = pmb->js;
    ib1 = pmb->is - 1;
    ib2 = pmb->ie + 1;
    kb1 = pmb->ks;
    kb2 = pmb->ke;
  }
  if (ox3 == 1) {
    DirNum = X3DIR;
    jb1 = pmb->js;
    jb2 = pmb->je;
    ib1 = pmb->is - 1;
    ib2 = pmb->ie + 1;
    kb1 = pmb->ke + 1;
    kb2 = pmb->ke + 1;
  }
  if (ox3 == -1) {
    DirNum = X3DIR;
    jb1 = pmb->js;
    jb2 = pmb->je;
    ib1 = pmb->is - 1;
    ib2 = pmb->ie + 1;
    kb1 = pmb->ks;
    kb2 = pmb->ks;
  }
  int dsize = ((kb2 - kb1 + 1) * (jb2 - jb1 + 1) * (ib2 - ib1 + 1) * NWAVE);
  Real *data = new Real[dsize];
  // Calculate the tag for receiving
  if (ox2 == -1)
    DirTag = 0 + 4 * tg_gid + 24 * (1 << (loc.level - 2)) * pmb->gid;
  if (ox2 == 1)
    DirTag = 1 + 4 * tg_gid + 24 * (1 << (loc.level - 2)) * pmb->gid;
  if (ox3 == -1)
    DirTag = 2 + 4 * tg_gid + 24 * (1 << (loc.level - 2)) * pmb->gid;
  if (ox3 == 1)
    DirTag = 3 + 4 * tg_gid + 24 * (1 << (loc.level - 2)) * pmb->gid;

#ifdef MPI_PARALLEL
  MPI_Recv(data, dsize, MPI_DOUBLE, tg_rank, DirTag, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
#endif

  // =======
  // int test;
  // probe MPI communications.  This is a bit of black magic that seems to
  // promote communications to top of stack and gets them to complete more
  // quickly MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
  //             MPI_STATUS_IGNORE);
  // =======
  // MPI_Test(&recv_request_[ownTag], &test, MPI_STATUS_IGNORE);

  int offset = 0;
  for (int n = 0; n < NWAVE; n++)
    for (int k = kb1; k <= kb2; k++)
      for (int j = jb1; j <= jb2; j++)
        for (int i = ib1; i <= ib2; i++)
          if (Left)
            L3DValues[DirNum](n, k, j, i) = data[offset++];
          else
            R3DValues[DirNum](n, k, j, i) = data[offset++];
  delete[] data;
}
