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
//! \param a $x = \tan(\xi)$ coordinate
//! \param b $y = \tan(\eta)$ coordinate
inline void vel_zab_p1_to_p2(Real *v1, Real *v2, Real *v3, Real a, Real b) {
  vel_zab_to_zxy(v1, v2, v3, a, b);

  // fill in code here

  vel_zxy_to_zab(v1, v2, v3, a, b);
}

}  // namespace CubedSphereUtility

#endif  // SRC_EXO3_VELOCITY_ROTATION_HPP_
