/** @file material_point.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Monday May 31, 2021 14:25:25 PDT
 * @bug No known bugs.
 */

// C/C++ headers
#include <iostream>

// Athena++ headers
#include "material_point.hpp"

MaterialPoint::MaterialPoint()
    : next(nullptr),
      id(0),
      type(0),
      time(0),
      rho(0),
      x1(0.),
      x2(0.),
      x3(0.),
      v1(0.),
      v2(0.),
      v3(0.) {
#if NREAL_PARTICLE_DATA > 0
  for (int i = 0; i < NREAL_PARTICLE_DATA; ++i) rr[i] = 0.;
#endif

#if NINT_PARTICLE_DATA > 0
  for (int i = 0; i < NINT_PARTICLE_DATA; ++i) ii[i] = 0;
#endif
}

MaterialPoint::~MaterialPoint() {}

MaterialPoint::MaterialPoint(MaterialPoint const& other) {
  if (this == &other) return;
  *this = other;
}

MaterialPoint& MaterialPoint::operator=(MaterialPoint const& other) {
  next = other.next;
  id = other.id;
  type = other.type;
  time = other.time;
  rho = other.rho;
  x1 = other.x1;
  x2 = other.x2;
  x3 = other.x3;
  v1 = other.v1;
  v2 = other.v2;
  v3 = other.v3;
#if NREAL_PARTICLE_DATA > 0
  for (int i = 0; i < NREAL_PARTICLE_DATA; ++i) rr[i] = other.rr[i];
#endif

#if NINT_PARTICLE_DATA > 0
  for (int i = 0; i < NINT_PARTICLE_DATA; ++i) ii[i] = other.ii[i];
#endif
  return *this;
}

std::ostream& operator<<(std::ostream& os, MaterialPoint const& pt) {
  os << "density: " << pt.rho << " type: " << pt.type << std::endl
     << "x1: " << pt.x1 << " v1: " << pt.v1 << std::endl
     << "x2: " << pt.x2 << " v2: " << pt.v2 << std::endl
     << "x3: " << pt.x3 << " v3: " << pt.v3 << std::endl;
#if NREAL_PARTICLE_DATA > 0
  os << "other real data: ";
  for (int i = 0; i < NREAL_PARTICLE_DATA; ++i) os << pt.rr[i] << " ";
  os << std::endl;
#endif

#if NINT_PARTICLE_DATA > 0
  os << "other integer data: ";
  for (int i = 0; i < NREAL_PARTICLE_DATA; ++i) os << pt.ii[i] << " ";
  os << std::endl;
#endif
  return os;
}
