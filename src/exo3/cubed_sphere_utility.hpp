#ifndef SRC_EXO3_CUBED_SPHERE_UTILITY_HPP_
#define SRC_EXO3_CUBED_SPHERE_UTILITY_HPP_

// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>

namespace cs {

void PackData(const AthenaArray<Real> &src, Real *buf, int sn, int en, int si,
              int ei, int sj, int ej, int sk, int ek, int &offset, int ox1,
              int ox2, int ox3, LogicalLocation const &loc, int TypeFlag);

// Helper functions adapted from Paul
void VecTransABPFromRLL(Real X, Real Y, int blockID, Real U, Real V, Real *V2,
                        Real *V3);
void VecTransRLLFromABP(Real X, Real Y, int blockID, Real V2, Real V3, Real *U,
                        Real *V);
void RLLFromXYP(Real dX, Real dY, int nP, Real &lon, Real &lat);
void XYPFromRLL(Real lon, Real lat, Real &dX, Real &dY, int &nP);

template <typename A>
void CovariantToContravariant(A a, Real cth) {
  Real v = a[IVY];
  Real w = a[IVZ];
  Real sth2 = 1. - cth * cth;

  a[IVY] = v / sth2 - w * cth / sth2;
  a[IVZ] = -v * cth / sth2 + w / sth2;
}

template <typename A>
void ContravariantToCovariant(A a, Real cth) {
  Real v = a[IVY];
  Real w = a[IVZ];
  a[IVY] = v + w * cth;
  a[IVZ] = w + v * cth;
}

void get_latlon_on_sphere(Real *lat, Real *lon, MeshBlock const *pmb, int k,
                          int j, int i);

}  // namespace cs

#endif  // SRC_EXO3_CUBED_SPHERE_UTILITY_HPP_
