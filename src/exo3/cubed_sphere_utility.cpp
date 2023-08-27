// C/C++
#include <algorithm>
#include <cmath>
#include <iostream>
#include <ostream>
#include <sstream>

// athena
#include <athena/bvals/bvals.hpp>
#include <athena/coordinates/coordinates.hpp>

// exo3
#include "cubed_sphere.hpp"
#include "cubed_sphere_utility.hpp"

namespace CubedSphereUtility {

void PackDataR3(const AthenaArray<Real> &src, Real *buf, int sn, int en, int si,
                int ei, int sj, int ej, int sk, int ek, int &offset) {
  for (int n = sn; n <= en; n++) {
    for (int j = sj; j <= ej; j++) {
      for (int k = ek; k >= sk; k--) {
#pragma omp simd
        for (int i = si; i <= ei; i++) buf[offset++] = src(n, k, j, i);
      }
    }
  }
  return;
}

void PackDataR2(const AthenaArray<Real> &src, Real *buf, int sn, int en, int si,
                int ei, int sj, int ej, int sk, int ek, int &offset) {
  for (int n = sn; n <= en; n++) {
    for (int k = ek; k >= sk; k--) {
      for (int j = ej; j >= sj; j--) {
#pragma omp simd
        for (int i = si; i <= ei; i++) buf[offset++] = src(n, k, j, i);
      }
    }
  }
  return;
}

void PackDataR1(const AthenaArray<Real> &src, Real *buf, int sn, int en, int si,
                int ei, int sj, int ej, int sk, int ek, int &offset) {
  for (int n = sn; n <= en; n++) {
    for (int j = ej; j >= sj; j--) {
      for (int k = sk; k <= ek; k++) {
#pragma omp simd
        for (int i = si; i <= ei; i++) buf[offset++] = src(n, k, j, i);
      }
    }
  }
  return;
}

void PackDataR0(const AthenaArray<Real> &src, Real *buf, int sn, int en, int si,
                int ei, int sj, int ej, int sk, int ek, int &offset) {
  for (int n = sn; n <= en; n++) {
    for (int k = sk; k <= ek; k++) {
      for (int j = sj; j <= ej; j++) {
#pragma omp simd
        for (int i = si; i <= ei; i++) buf[offset++] = src(n, k, j, i);
      }
    }
  }
  return;
}

Real CalculateInterpLocations(int loc_n, int N_blk, int k, bool GhostZone) {
  // Calculate the angular locations of the ghost zone interpolations
  // N_blk: the total number of points along the boundary
  // k: levels into the ghost zone (0, 1, 2 in interior = -1, -2, -3 in ghost)
  // loc_n: local index in the panel (-N_blk/2 to N_blk/2+1)
  Real temp0 = tan((Real)1.0 * (N_blk - 1 - 2 * k) / (N_blk * 4) * PI);
  if (GhostZone) {
    temp0 = 1.0 / temp0;
    temp0 = 1.0 + temp0 * temp0;
  } else {
    temp0 = 1.0 + temp0 * temp0;
  }
  if (loc_n < 0)
    return -acos(sqrt(temp0 / (temp0 + pow(tan((Real)1.0 * (2 * loc_n + 1) /
                                               (N_blk * 4) * PI),
                                           2.0))));
  else
    return acos(sqrt(temp0 / (temp0 + pow(tan((Real)1.0 * (2 * loc_n + 1) /
                                              (N_blk * 4) * PI),
                                          2.0))));
}

void InteprolateX2(const AthenaArray<Real> &src, AthenaArray<Real> &tgt,
                   LogicalLocation const &loc, int DirInv, int TgtSide,
                   int TgtID, int sn, int en, int si, int ei, int sj, int ej,
                   int sk, int ek) {
  // Interpolation along X2 (j) axis, used before sending data to X3 (k) axis
  // Get the local indices
  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int local_lx2 = loc.lx2 - (lv2_lx2 << (loc.level - 2));
  int bound_lim = (1 << (loc.level - 2)) - 1;
  int meshblock_size = ej - sj + 1;
  int N_blk = meshblock_size *
              (bound_lim + 1);  // N in X2 direction for each panel. This value
                                // is required to be an even number
  int n_start = local_lx2 * meshblock_size - N_blk / 2;
  Real *src_x2 = new Real[N_blk];  // Need to calculate source indices along the
                                   // whole panel boundary
  Real *src_coord = new Real[N_blk];  // The src coordinate along panel boundary
  Real *tgt_coord =
      new Real[ej - sj + 1];  // The tgt coordinate along panel boundary
  Real *tgt_x2 = new Real[ej - sj + 1];
  int SrcSide;

  for (int n = sn; n <= en; n++) {
    for (int k = sk; k <= ek; k++) {
      int k_now;
      if (sk < ej - NGHOST) {  // Calculate the location into the ghost zones
        k_now = k - sk;
        SrcSide = -1;
      } else {
        k_now = ek - k;
        SrcSide = 1;
      }

#pragma omp simd
      for (int i = si; i <= ei; i++) {
        for (int j = 0; j < N_blk; j++) {
          // Calculate coefficients for src, need to go from 0 to N_blk to cover
          // all
          src_x2[j] =
              CalculateInterpLocations(j - N_blk / 2, N_blk, k_now, false);
          // Calculate the coordinate locations, going from -pi/4 to pi/4
          src_coord[j] = PI / 2.0 / N_blk * (j + 0.5) - PI / 4.0;
        }
        for (int j = sj; j <= ej; j++) {
          // Calculate coefficients for tgt
          tgt_x2[j - sj] =
              CalculateInterpLocations(n_start + j - sj, N_blk, k_now, true);
          if (DirInv == 1) {
            tgt_coord[j - sj] = -src_coord[n_start + N_blk / 2 + j - sj];
          } else {
            tgt_coord[j - sj] = src_coord[n_start + N_blk / 2 + j - sj];
          }
        }
        int src_pointer = 0;
        for (int j = sj; j <= ej; j++) {
          // Interpolate to target array, linear interpolation used here
          while (tgt_x2[j - sj] > src_x2[src_pointer + 1])
            src_pointer++;  // Find the left location of src

          Real y1 = src_x2[src_pointer];
          Real y2 = src_x2[src_pointer + 1];
          Real yq = tgt_x2[j - sj];
          if (n == IVY || n == IVZ) {
            // Projection needed, find the tgt locations first
            Real v1y =
                src(IVY, k, src_pointer - n_start - N_blk / 2 + sj,
                    i);  // When pointing at n_start+N_blk, we are looking at sj
            Real v1z = src(IVZ, k, src_pointer - n_start - N_blk / 2 + sj, i);
            Real v2y =
                src(IVY, k, src_pointer + 1 - n_start - N_blk / 2 + sj, i);
            Real v2z =
                src(IVZ, k, src_pointer + 1 - n_start - N_blk / 2 + sj, i);
            Real tgt_cy = tan(tgt_coord[j - sj]);
            Real tgt_cz =
                tan(TgtSide * (PI / 2.0 / (N_blk) * (k_now + 0.5) + PI / 4.0));
            Real src_cy1 = src_coord[src_pointer];
            Real src_cy2 = src_coord[src_pointer + 1];
            Real src_cz =
                tan(SrcSide * (-PI / 2.0 / (N_blk) * (k_now + 0.5) + PI / 4.0));
            Real vy = ((y2 - yq) * v1y + (yq - y1) * v2y) / (y2 - y1);
            Real vz = ((y2 - yq) * v1z + (yq - y1) * v2z) / (y2 - y1);
            Real src_cy =
                tan(((y2 - yq) * src_cy1 + (yq - y1) * src_cy2) / (y2 - y1));
            int blockID = CubedSphere::FindBlockID(loc);
            Real U, V;
            // Raise the vectors to contravariant
            Real sth = sqrt(1 + src_cy * src_cy + src_cz * src_cz) /
                       sqrt(1 + src_cy * src_cy) / sqrt(1 + src_cz * src_cz);
            Real cth = -src_cy * src_cz / sqrt(1 + src_cy * src_cy) /
                       sqrt(1 + src_cz * src_cz);
            Real vyt = vy;
            Real vzt = vz;
            vy = vyt / sth / sth - vzt * cth / sth / sth;
            vz = vzt / sth / sth - vyt * cth / sth / sth;

            // Trasform to global coordinate and back
            VecTransRLLFromABP(src_cy, src_cz, blockID, vy, vz, &U, &V);
            VecTransABPFromRLL(tgt_cy, tgt_cz, TgtID, U, V, &vy, &vz);

            // Lower the vectors to covariant
            sth = sqrt(1 + tgt_cy * tgt_cy + tgt_cz * tgt_cz) /
                  sqrt(1 + tgt_cy * tgt_cy) / sqrt(1 + tgt_cz * tgt_cz);
            cth = -tgt_cy * tgt_cz / sqrt(1 + tgt_cy * tgt_cy) /
                  sqrt(1 + tgt_cz * tgt_cz);
            vyt = vy;
            vzt = vz;
            vy = vyt + vzt * cth;
            vz = vzt + vyt * cth;

            if (n == IVY) {
              tgt(n - sn, k - sk, j - sj, i - si) = vy;
            } else {  // n==IVZ
              tgt(n - sn, k - sk, j - sj, i - si) = vz;
            }
          } else {
            Real v1 = src(n, k, src_pointer - n_start - N_blk / 2 + sj, i);
            Real v2 = src(n, k, src_pointer + 1 - n_start - N_blk / 2 + sj, i);
            tgt(n - sn, k - sk, j - sj, i - si) =
                ((y2 - yq) * v1 + (yq - y1) * v2) / (y2 - y1);
          }
        }
      }
    }
  }
  delete[] src_x2;
  delete[] tgt_x2;
  delete[] src_coord;
  delete[] tgt_coord;  // Release memory
}

void InteprolateX3(const AthenaArray<Real> &src, AthenaArray<Real> &tgt,
                   LogicalLocation const &loc, int DirInv, int TgtSide,
                   int TgtID, int sn, int en, int si, int ei, int sj, int ej,
                   int sk, int ek) {
  // Interpolation along X3 (k) axis, used before sending data to ghost zone in
  // X2 (j) direction Get the local indices
  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int local_lx3 = loc.lx3 - (lv2_lx3 << (loc.level - 2));
  int bound_lim = (1 << (loc.level - 2)) - 1;
  int meshblock_size = ek - sk + 1;
  int N_blk = meshblock_size *
              (bound_lim + 1);  // N in X2 direction for each panel. This value
                                // is required to be an even number
  int n_start = local_lx3 * meshblock_size - N_blk / 2;
  Real *src_x3 = new Real[N_blk];  // Need to calculate source indices along the
                                   // whole panel boundary
  Real *src_coord = new Real[N_blk];  // The src coordinate along panel boundary
  Real *tgt_coord =
      new Real[ek - sk + 1];  // The tgt coordinate along panel boundary
  Real *tgt_x3 = new Real[ek - sk + 1];
  int SrcSide;

  for (int n = sn; n <= en; n++) {
    for (int j = sj; j <= ej; j++) {
      int k_now;
      if (sj < ek - NGHOST) {  // Calculate the location into the ghost zones
        k_now = j - sj;
        SrcSide = -1;
      } else {
        k_now = ej - j;
        SrcSide = 1;
      }
#pragma omp simd
      for (int i = si; i <= ei; i++) {
        for (int k = 0; k < N_blk; k++) {
          // Calculate coefficients for src, need to go from 0 to N_blk to cover
          // all
          src_x3[k] =
              CalculateInterpLocations(k - N_blk / 2, N_blk, k_now, false);
          // Calculate the coordinate locations, going from -pi/4 to pi/4
          src_coord[k] = PI / 2.0 / N_blk * (k + 0.5) - PI / 4.0;
        }
        for (int k = sk; k <= ek; k++) {
          // Calculate coefficients for tgt
          tgt_x3[k - sk] =
              CalculateInterpLocations(n_start + k - sk, N_blk, k_now, true);
          if (DirInv == 1) {
            tgt_coord[k - sk] = -src_coord[n_start + N_blk / 2 + k - sk];
          } else {
            tgt_coord[k - sk] = src_coord[n_start + N_blk / 2 + k - sk];
          }
        }
        int src_pointer = 0;
        for (int k = sk; k <= ek; k++) {
          // Interpolate to target array, linear interpolation used here
          while (tgt_x3[k - sk] > src_x3[src_pointer + 1])
            src_pointer++;  // Find the left location of src
          Real y1 = src_x3[src_pointer];
          Real y2 = src_x3[src_pointer + 1];
          Real yq = tgt_x3[k - sk];
          if (n == IVY || n == IVZ) {
            // Projection needed, find the tgt locations first
            Real v1y = src(IVY, src_pointer - n_start - N_blk / 2 + sk, j, i);
            Real v1z = src(IVZ, src_pointer - n_start - N_blk / 2 + sk, j, i);
            Real v2y =
                src(IVY, src_pointer + 1 - n_start - N_blk / 2 + sk, j, i);
            Real v2z =
                src(IVZ, src_pointer + 1 - n_start - N_blk / 2 + sk, j, i);
            Real tgt_cz = tan(tgt_coord[k - sk]);
            Real tgt_cy =
                tan(TgtSide * (PI / 2.0 / (N_blk) * (k_now + 0.5) + PI / 4.0));
            Real src_cz1 = src_coord[src_pointer];
            Real src_cz2 = src_coord[src_pointer + 1];
            Real src_cy =
                tan(SrcSide * (-PI / 2.0 / (N_blk) * (k_now + 0.5) + PI / 4.0));
            Real vy = ((y2 - yq) * v1y + (yq - y1) * v2y) / (y2 - y1);
            Real vz = ((y2 - yq) * v1z + (yq - y1) * v2z) / (y2 - y1);
            Real src_cz =
                tan(((y2 - yq) * src_cz1 + (yq - y1) * src_cz2) / (y2 - y1));
            int blockID = CubedSphere::FindBlockID(loc);
            Real U, V;

            // Raise the vectors to contravariant
            Real sth = sqrt(1 + src_cy * src_cy + src_cz * src_cz) /
                       sqrt(1 + src_cy * src_cy) / sqrt(1 + src_cz * src_cz);
            Real cth = -src_cy * src_cz / sqrt(1 + src_cy * src_cy) /
                       sqrt(1 + src_cz * src_cz);
            Real vyt = vy;
            Real vzt = vz;
            vy = vyt / sth / sth - vzt * cth / sth / sth;
            vz = vzt / sth / sth - vyt * cth / sth / sth;

            // Trasform to global coordinate and back
            VecTransRLLFromABP(src_cy, src_cz, blockID, vy, vz, &U, &V);
            if (((TgtID == 1) && (blockID == 3)) ||
                ((TgtID == 1) && (blockID == 4)) ||
                ((TgtID == 5) && (blockID == 3)) ||
                ((TgtID == 5) && (blockID == 4)) ||
                ((TgtID == 3) && (blockID == 1)) ||
                ((TgtID == 4) && (blockID == 1)) ||
                ((TgtID == 3) && (blockID == 5)) ||
                ((TgtID == 4) && (blockID == 5))) {
              Real tmp = tgt_cy;
              tgt_cy = tgt_cz;
              tgt_cz = tmp;
            }
            VecTransABPFromRLL(tgt_cy, tgt_cz, TgtID, U, V, &vy, &vz);

            // Lower the vectors to covariant
            sth = sqrt(1 + tgt_cy * tgt_cy + tgt_cz * tgt_cz) /
                  sqrt(1 + tgt_cy * tgt_cy) / sqrt(1 + tgt_cz * tgt_cz);
            cth = -tgt_cy * tgt_cz / sqrt(1 + tgt_cy * tgt_cy) /
                  sqrt(1 + tgt_cz * tgt_cz);
            vyt = vy;
            vzt = vz;
            vy = vyt + vzt * cth;
            vz = vzt + vyt * cth;

            if (n == IVY) {
              tgt(n - sn, k - sk, j - sj, i - si) = vy;
            } else {  // n==IVZ
              tgt(n - sn, k - sk, j - sj, i - si) = vz;
            }
            if (n == IVY) {
              tgt(n - sn, k - sk, j - sj, i - si) = vy;
            } else {  // n==IVZ
              tgt(n - sn, k - sk, j - sj, i - si) = vz;
            }
          } else {
            if ((src_pointer - n_start - N_blk / 2 + sk > ek) ||
                (src_pointer - n_start - N_blk / 2 + sk < sk)) {
            }
            Real v1 = src(n, src_pointer - n_start - N_blk / 2 + sk, j, i);
            Real v2 = src(n, src_pointer + 1 - n_start - N_blk / 2 + sk, j, i);
            tgt(n - sn, k - sk, j - sj, i - si) =
                ((y2 - yq) * v1 + (yq - y1) * v2) / (y2 - y1);
          }
        }
      }
    }
  }
  delete[] src_x3;
  delete[] tgt_x3;
  delete[] src_coord;
  delete[] tgt_coord;  // Release memory
}

void PackData(const AthenaArray<Real> &src, Real *buf, int sn, int en, int si,
              int ei, int sj, int ej, int sk, int ek, int &offset, int ox1,
              int ox2, int ox3, LogicalLocation const &loc) {
  // Find the block ID
  int blockID = CubedSphere::FindBlockID(loc);

  // Table of wether direction inversion is needed
  const int dinv[6][4] = {// To access: dinv[source_id][target_dir]
                          {1, 0, 0, 1}, {0, 0, 0, 0}, {0, 1, 0, 0},
                          {1, 0, 0, 0}, {0, 1, 1, 0}, {1, 1, 0, 0}};

  // Table of which side (+pi/4 or -pi/4) is in touch with this panel
  const int tgside[6][4] = {// To access: dinv[source_id][target_dir]
                            {-1, -1, -1, -1}, {1, -1, 1, -1}, {-1, -1, 1, -1},
                            {1, 1, 1, -1},    {1, 1, 1, 1},   {-1, 1, 1, -1}};

  // Table of which panel is in touch with this panel
  const int tgbid[6][4] = {// To access: tgbid[source_id][target_dir]
                           {6, 2, 3, 4}, {1, 5, 3, 4}, {1, 5, 6, 2},
                           {1, 5, 2, 6}, {2, 6, 3, 4}, {1, 5, 4, 3}};

  // Bypass Corner cases
  if ((ox2 + ox3 == 0) || (ox2 + ox3 == 2) || (ox2 + ox3 == -2) || (ox1 != 0)) {
    PackDataR0(src, buf, sn, en, si, ei, sj, ej, sk, ek, offset);
    return;
  }

  // Get the local indices
  int lv2_lx2 = loc.lx2 >> (loc.level - 2);
  int lv2_lx3 = loc.lx3 >> (loc.level - 2);
  int local_lx2 = loc.lx2 - (lv2_lx2 << (loc.level - 2));
  int local_lx3 = loc.lx3 - (lv2_lx3 << (loc.level - 2));
  int bound_lim = (1 << (loc.level - 2)) - 1;

  // Work on interpolation
  AthenaArray<Real> interpolatedSrc;
  interpolatedSrc.NewAthenaArray(en - sn + 1, ek - sk + 1, ej - sj + 1,
                                 ei - si + 1);

  switch (blockID) {
    case 1:
      if (ox3 == 1 && local_lx3 == bound_lim) {
        InteprolateX2(src, interpolatedSrc, loc, dinv[blockID - 1][3],
                      tgside[blockID - 1][3], tgbid[blockID - 1][3], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR1(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      if (ox2 == 1 && local_lx2 == bound_lim) {
        InteprolateX3(src, interpolatedSrc, loc, dinv[blockID - 1][1],
                      tgside[blockID - 1][1], tgbid[blockID - 1][1], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR0(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      if (ox2 == -1 && local_lx2 == 0) {
        InteprolateX3(src, interpolatedSrc, loc, dinv[blockID - 1][0],
                      tgside[blockID - 1][0], tgbid[blockID - 1][0], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR2(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      if (ox3 == -1 && local_lx3 == 0) {
        InteprolateX2(src, interpolatedSrc, loc, dinv[blockID - 1][2],
                      tgside[blockID - 1][2], tgbid[blockID - 1][2], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR3(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      break;
    case 2:
      if (ox3 == 1 && local_lx3 == bound_lim) {
        InteprolateX2(src, interpolatedSrc, loc, dinv[blockID - 1][3],
                      tgside[blockID - 1][3], tgbid[blockID - 1][3], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR0(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      if (ox2 == 1 && local_lx2 == bound_lim) {
        InteprolateX3(src, interpolatedSrc, loc, dinv[blockID - 1][1],
                      tgside[blockID - 1][1], tgbid[blockID - 1][1], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR0(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      if (ox2 == -1 && local_lx2 == 0) {
        InteprolateX3(src, interpolatedSrc, loc, dinv[blockID - 1][0],
                      tgside[blockID - 1][0], tgbid[blockID - 1][0], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR0(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      if (ox3 == -1 && local_lx3 == 0) {
        InteprolateX2(src, interpolatedSrc, loc, dinv[blockID - 1][2],
                      tgside[blockID - 1][2], tgbid[blockID - 1][2], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR0(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      break;
    case 3:
      if (ox2 == 1 && local_lx2 == bound_lim) {
        InteprolateX3(src, interpolatedSrc, loc, dinv[blockID - 1][1],
                      tgside[blockID - 1][1], tgbid[blockID - 1][1], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR3(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      if (ox2 == -1 && local_lx2 == 0) {
        InteprolateX3(src, interpolatedSrc, loc, dinv[blockID - 1][0],
                      tgside[blockID - 1][0], tgbid[blockID - 1][0], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR1(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      if (ox3 == 1 && local_lx3 == bound_lim) {
        InteprolateX2(src, interpolatedSrc, loc, dinv[blockID - 1][3],
                      tgside[blockID - 1][3], tgbid[blockID - 1][3], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR0(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      if (ox3 == -1 && local_lx3 == 0) {
        InteprolateX2(src, interpolatedSrc, loc, dinv[blockID - 1][2],
                      tgside[blockID - 1][2], tgbid[blockID - 1][2], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR0(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      break;
    case 4:
      if (ox2 == 1 && local_lx2 == bound_lim) {
        InteprolateX3(src, interpolatedSrc, loc, dinv[blockID - 1][1],
                      tgside[blockID - 1][1], tgbid[blockID - 1][1], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR1(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      if (ox2 == -1 && local_lx2 == 0) {
        InteprolateX3(src, interpolatedSrc, loc, dinv[blockID - 1][0],
                      tgside[blockID - 1][0], tgbid[blockID - 1][0], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR3(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      if (ox3 == 1 && local_lx3 == bound_lim) {
        InteprolateX2(src, interpolatedSrc, loc, dinv[blockID - 1][3],
                      tgside[blockID - 1][3], tgbid[blockID - 1][3], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR0(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      if (ox3 == -1 && local_lx3 == 0) {
        InteprolateX2(src, interpolatedSrc, loc, dinv[blockID - 1][2],
                      tgside[blockID - 1][2], tgbid[blockID - 1][2], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR0(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      break;
    case 5:
      if (ox3 == 1 && local_lx3 == bound_lim) {
        InteprolateX2(src, interpolatedSrc, loc, dinv[blockID - 1][3],
                      tgside[blockID - 1][3], tgbid[blockID - 1][3], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR3(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      if (ox2 == 1 && local_lx2 == bound_lim) {
        InteprolateX3(src, interpolatedSrc, loc, dinv[blockID - 1][1],
                      tgside[blockID - 1][1], tgbid[blockID - 1][1], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR2(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      if (ox3 == -1 && local_lx3 == 0) {
        InteprolateX2(src, interpolatedSrc, loc, dinv[blockID - 1][2],
                      tgside[blockID - 1][2], tgbid[blockID - 1][2], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR1(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      if (ox2 == -1 && local_lx2 == 0) {
        InteprolateX3(src, interpolatedSrc, loc, dinv[blockID - 1][0],
                      tgside[blockID - 1][0], tgbid[blockID - 1][0], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR0(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      break;
    case 6:
      if (ox2 == 1 && local_lx2 == bound_lim) {
        InteprolateX3(src, interpolatedSrc, loc, dinv[blockID - 1][1],
                      tgside[blockID - 1][1], tgbid[blockID - 1][1], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR2(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      if (ox2 == -1 && local_lx2 == 0) {
        InteprolateX3(src, interpolatedSrc, loc, dinv[blockID - 1][0],
                      tgside[blockID - 1][0], tgbid[blockID - 1][0], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR2(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      if (ox3 == 1 && local_lx3 == bound_lim) {
        InteprolateX2(src, interpolatedSrc, loc, dinv[blockID - 1][3],
                      tgside[blockID - 1][3], tgbid[blockID - 1][3], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR0(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
      if (ox3 == -1 && local_lx3 == 0) {
        InteprolateX2(src, interpolatedSrc, loc, dinv[blockID - 1][2],
                      tgside[blockID - 1][2], tgbid[blockID - 1][2], sn, en, si,
                      ei, sj, ej, sk, ek);
        PackDataR0(interpolatedSrc, buf, 0, en - sn, 0, ei - si, 0, ej - sj, 0,
                   ek - sk, offset);
        return;
      }
  }
  PackDataR0(src, buf, sn, en, si, ei, sj, ej, sk, ek, offset);
  return;
}

}  // namespace CubedSphereUtility
