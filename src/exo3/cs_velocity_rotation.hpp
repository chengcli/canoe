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

}

//! Transform local cartesian velocity to cubed sphere velocity
inline void vel_zxy_to_zab(Real *v1, Real *v2, Real *v3, Real a, Real b) {
    Real vx = v2;
    Real vy = v3;
    Real delta = pow(pow(vx,2)+pow(vy,2)+1,1/2);
    Real C = pow(1+pow(x,2),1/2);
    Real D = pow(1+pow(y,2),1/2);
    v1 = (vz+pow(vx,2)+pow(vy,2))/delta;
    v2 = (-vx*vz/D+vx*(1+pow(vy,2))/D-vx*pow(vy,2)/D)/delta;
    v3 = (-vy*vz/C-pow(vx,2)*vy/C+(1+pow(vx,2))*vy/C)/delta;
}

//! Transform cubed sphere velocity from panel 1 to panel 2
//! \param a $x = \tan(\xi)$ coordinates
//! \param b $y = \tan(\eta)$ coordinat
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

inline void vel_zab_from_p1_test(std::vector<Real> &v, int panel) {
    Real vz = v[0];
    Real vx = v[1];
    Real vy = v[2];
    switch (panel) {
    case 2:
      // z->y, x->-x, y->z
      //(*vx) *= -1;
      v = {vy, -vx, vz};
      break;
    case 3:
      // z->-x, x->z, y->y
      //(*vz) *= -1;
      v = {vx, -vz, vy};
      break;
    case 4:
      // z->-x, x->-y, y->z
      //(*vx) *= -1;
      //(*vy) *= -1;
      v = {-vx, -vy, vz};
      break;
    case 6:
      // z->-y, x->x, y->z
      //(*vy) *= -1;
      v = {-vy, vx, vz};
      break;
  }
 }
}  // namespace CubedSphereUtility

#endif  // SRC_EXO3_VELOCITY_ROTATION_HPP_
