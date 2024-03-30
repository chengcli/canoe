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
inline void vel_zab_to_zxy(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  Real x = tan(a);
  Real y = tan(b);

  Real vx = *v2;
  Real vy = *v3;
  Real vz = *v1;

  Real delta = sqrt(x*x + y*y + 1);
  Real C = sqrt(1 + x*x);
  Real D = sqrt(1 + y*y);

  *v1 = (vz - D * x * vx - C * y * vy) / delta;
  *v2 = 
      (x * vz + D * vx) / delta;
  *v3 =
      (y * vz + C * vy) / delta;
}

//! Transform local cartesian velocity to cubed sphere velocity
inline void vel_zxy_to_zab(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  Real x = tan(a);
  Real y = tan(b);

  Real vx = *v2;
  Real vy = *v3;
  Real vz = *v1;

  Real delta = sqrt(x*x + y*y + 1);
  Real C = sqrt(1 + x*x);
  Real D = sqrt(1 + y*y);

  *v1 = (vz + x * vx + y * vy) / delta;
  *v2 =
      (-x * vz / D + vx * (1 + y*y) / D - vy * x * y / D) / delta;
  *v3 =
      (-y * vz / C - x * y * vx / C + (1 + x*x) * vy / C) / delta;
}

//! Transform cubed sphere velocity from panel 1 to panel 2
//! \param a $x = \tan(\xi)$ coordinates
//! \param b $y = \tan(\eta)$ coordinat
inline void vel_zab_from_p1(Real *vz, Real *vx, Real *vy, Real a, Real b,
                            int panel) {
  vel_zab_to_zxy(vz, vx, vy, a, b);
  Real v1 = *vz;
  Real v2 = *vx;
  Real v3 = *vy;

  switch (panel) {
    case 2:
      // z->y, x->-x, y->z
      //(*vx) *= -1;
      //vel_zxy_to_zab(vy, vx, vz, a, b);
      *vz = v3;
      *vx = -v2;
      *vy = v1;
      vel_zxy_to_zab(vz, vx, vy, a, b);
      break;
    case 3:
      // z->-x, x->z, y->y
      //(*vz) *= -1;
      //vel_zxy_to_zab(vx, vz, vy, a, b);
      *vz = v2;
      *vx = -v1;
      *vy = v3;
      vel_zxy_to_zab(vz, vx, vy, a, b);
      break;
    case 4:
      // z->-x, x->-y, y->z
      //(*vx) *= -1;
      //(*vy) *= -1;
      //vel_zxy_to_zab(vx, vy, vz, a, b);
      *vz = -v2;
      *vx = -v3;
      *vy = v1;
      vel_zxy_to_zab(vz, vx, vy, a, b);
      break;
    case 6:
      // z->-y, x->x, y->z
      //(*vy) *= -1;
      //vel_zxy_to_zab(vy, vx, vz, a, b);
      *vz = -v3;
      *vx = v2;
      *vy = v1;
      vel_zxy_to_zab(vz, vx, vy, a, b);
      break;
  }
}

inline void vel_zab_from_p2(Real *vz, Real *vx, Real *vy, Real a, Real b,
                            int panel) {
  vel_zab_to_zxy(vz, vx, vy, a, b);
  Real v1 = *vz;
  Real v2 = *vx;
  Real v3 = *vy;
  switch (panel) {
    case 1:
      // z->y, x->-x, y->z
      *vz = v3;
      *vx = -v2;
      *vy = v1;
      vel_zxy_to_zab(vz, vx, vy, a, b);
      break;
    case 2:
      break;
    case 4:
      break;
    case 5:
      break;
  }
}

inline void vel_zab_from_p3(Real *vz, Real *vx, Real *vy, Real a, Real b,
                            int panel) {
  vel_zab_to_zxy(vz, vx, vy, a, b);
  Real v1 = *vz;
  Real v2 = *vx;
  Real v3 = *vy;
  switch (panel) {
    case 1:
      // z->x, x->-z, y->y
      *vz = -v2;
      *vx = v1;
      *vy = v3;
      vel_zxy_to_zab(vz, vx, vy, a, b);
      break;
    case 2:
      break;
    case 5:
      break;
    case 6:
      break;
  }
}

inline void vel_zab_from_p4(Real *vz, Real *vx, Real *vy, Real a, Real b,
                            int panel) {
  vel_zab_to_zxy(vz, vx, vy, a, b);
  Real v1 = *vz;
  Real v2 = *vx;
  Real v3 = *vy;
  switch (panel) {
    case 1:
      // z->y, x->-z, y->-x
      *vz = v3;
      *vx = -v1;
      *vy = -v2;
      vel_zxy_to_zab(vz, vx, vy, a, b);
      break;
    case 2:
      break;
    case 5:
      break;
    case 6:
      break;
  }
}

inline void vel_zab_from_p5(Real *vz, Real *vx, Real *vy, Real a, Real b,
                            int panel) {
  vel_zab_to_zxy(vz, vx, vy, a, b);
  Real v1 = *vz;
  Real v2 = *vx;
  Real v3 = *vy;
  switch (panel) {
    case 2:
      break;
    case 3:
      break;
    case 4:
      break;
    case 6:
      break;
  }
}

inline void vel_zab_from_p6(Real *vz, Real *vx, Real *vy, Real a, Real b,
                            int panel) {
  vel_zab_to_zxy(vz, vx, vy, a, b);
  Real v1 = *vz;
  Real v2 = *vx;
  Real v3 = *vy;
  switch (panel) {
    case 1:
      // z->y, x->x, y->-z
      *vz = v3;
      *vx = v2;
      *vy = -v1;
      vel_zxy_to_zab(vz, vx, vy, a, b);
      break;
    case 3:
      break;
    case 4:
      break;
    case 5:
      break;
  }
}
}  // namespace CubedSphereUtility

#endif  // SRC_EXO3_VELOCITY_ROTATION_HPP_