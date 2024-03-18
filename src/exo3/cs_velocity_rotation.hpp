#ifndef SRC_EXO3_VELOCITY_ROTATION_HPP_
#define SRC_EXO3_VELOCITY_ROTATION_HPP_

#include <athena/athena.hpp>  // Real

//! Panel number assiangments
//! (1) Top
//! (2) Front
//! (3) Left
//!
//!                                           z
//!                ___________                |
//!                |\        .\           x___| (1)
//!                | \   1   . \              /
//!                |  \_________\           y/
//!                | 3 |     .  |
//!                \. .|......  |
//! z___ (3)        \  |    2 . |
//!    /|            \ |       .|             y
//!   / |             \|________|             |
//!  y  x                                 (2) |___x
//!                                          /
//!                                        z/

//! Panel number assiangments
//! (4) Right
//! (5) Bottom
//! (6) Back
//!                                        y  x
//!                __________              | /
//!                |\       .\             |/___z
//!                | \      . \           (4)
//!     y  z       |  \________\
//!     | /        |  |  6  .  |
//! x___|/         |..|......  |
//!      (6)       \  |     . 4|       (5) ___ x
//!                 \ |  5   . |          /|
//!                  \|_______.|         / |
//!                                     y  z

namespace CubedSphereUtility {

//! Transform cubed sphere velocity to local cartesian velocity
inline void vel_zab_to_zxy(Real *v1, Real *v2, Real *v3, Real a, Real b) {}

//! Transform local cartesian velocity to cubed sphere velocity
inline void vel_zxy_to_zab(Real *v1, Real *v2, Real *v3, Real a, Real b) {}

//! Transform cubed sphere velocity from panel 1 to panel 2
//! \param a $x = \tan(\xi)$ coordinates
//! \param b $y = \tan(\eta)$ coordinate
inline void vel_zab_from_p1(Real *vz, Real *vx, Real *vy, Real a, Real b,
                            int panel) {
  vel_zab_to_zxy(vz, vx, vy, a, b);

  switch (panel) {
    case 2:
      // z->y, x->-x, y->z
      (*vx) *= -1;
      vel_zxy_to_zab(vy, vx, vz, a, b);
      break;
    case 3:
      // z->-x, x->z, y->y
      (*vz) *= -1;
      vel_zxy_to_zab(vx, vz, vy, a, b);
      break;
    case 4:
      // z->-x, x->-y, y->z
      (*vx) *= -1;
      (*vy) *= -1;
      vel_zxy_to_zab(vx, vy, vz, a, b);
      break;
    case 6:
      // z->-y, x->x, y->z
      (*vy) *= -1;
      vel_zxy_to_zab(vy, vx, vz, a, b);
      break;
  }
}

inline void vel_zab_p2_to_p4(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {0, 0, 1}, {0, 1, 0}, {-1, 0, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p3_to_p5(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {0, 0, 1}, {0, 1, 0}, {-1, 0, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p4_to_p6(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {0, 0, 1}, {0, 1, 0}, {-1, 0, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p3_to_p1(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {0, 0, -1}, {0, 1, 0}, {1, 0, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p5_to_p3(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {0, 0, -1}, {0, 1, 0}, {1, 0, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p4_to_p2(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {0, 0, -1}, {0, 1, 0}, {1, 0, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p6_to_p4(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {0, 0, -1}, {0, 1, 0}, {1, 0, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p1_to_p5(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {-1, 0, 0}, {0, 1, 0}, {0, 0, -1}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p2_to_p6(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {-1, 0, 0}, {0, 1, 0}, {0, 0, -1}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p6_to_p2(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {-1, 0, 0}, {0, 1, 0}, {0, 0, -1}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p5_to_p1(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {-1, 0, 0}, {0, 1, 0}, {0, 0, -1}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p5_to_p2(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {-1, 0, 0}, {0, 0, -1}, {0, -1, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p2_to_p5(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {-1, 0, 0}, {0, -1, 0}, {0, 0, -1}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p5_to_p4(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {0, -1, 0}, {0, 0, -1}, {1, 0, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p4_to_p5(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {0, 0, 1}, {-1, 0, 0}, {0, -1, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p5_to_p6(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {1, 0, 0}, {0, 0, -1}, {0, 1, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p6_to_p5(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {1, 0, 0}, {0, 0, 1}, {0, -1, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p3_to_p2(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {0, 0, -1}, {1, 0, 0}, {0, -1, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p2_to_p3(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {0, 1, 0}, {0, 0, -1}, {-1, 0, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p3_to_p4(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {0, -1, 0}, {1, 0, 0}, {0, 0, 1}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p4_to_p3(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {0, 1, 0}, {-1, 0, 0}, {0, 0, 1}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p3_to_p6(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {0, 0, 1}, {1, 0, 0}, {0, 1, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p6_to_p3(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {0, 1, 0}, {0, 0, 1}, {1, 0, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p1_to_p2(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {0, 0, 1}, {1, 0, 0}, {0, 1, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p2_to_p1(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {0, 1, 0}, {0, 0, 1}, {1, 0, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p1_to_p4(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {0, -1, 0}, {0, 0, 1}, {-1, 0, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p4_to_p1(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {0, 0, -1}, {-1, 0, 0}, {0, 1, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p1_to_p6(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {-1, 0, 0}, {0, 0, 1}, {0, 1, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

inline void vel_zab_p6_to_p1(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
  M[3][3] = {
      {-1, 0, 0}, {0, 0, 1}, {0, 1, 0}};  // define the transformation matrix
  // multiply the velocity vector on the right of the transformation
  v1 = M[0][0] * v1 + M[1][0] * v2 + M[2][0] * v3;
  v2 = M[0][1] * v1 + M[1][1] * v2 + M[2][1] * v3;
  v2 = M[0][2] * v1 + M[1][2] * v2 + M[2][2] * v3;
  vel_zxy_to_zab(
      v1, v2, v3, a,
      b);  // translate the cubed sphere coordinate into cartesian coordinate
}

}  // namespace CubedSphereUtility

#endif  // SRC_EXO3_VELOCITY_ROTATION_HPP_
