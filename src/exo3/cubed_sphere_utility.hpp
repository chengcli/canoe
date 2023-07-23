#ifndef SRC_EXO3_CUBED_SPHERE_UTILITY_HPP_
#define SRC_EXO3_CUBED_SPHERE_UTILITY_HPP_

// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>

namespace CubedSphereUtility {

void PackData(const AthenaArray<Real> &src, Real *buf, int sn, int en, int si,
              int ei, int sj, int ej, int sk, int ek, int &offset, int ox1,
              int ox2, int ox3, LogicalLocation const &loc);

// Helper functions adapted from Paul
void VecTransABPFromRLL(Real X, Real Y, int blockID, Real U, Real V, Real *V2,
                        Real *V3);
void VecTransRLLFromABP(Real X, Real Y, int blockID, Real V2, Real V3, Real *U,
                        Real *V);
void RLLFromXYP(Real dX, Real dY, int nP, Real &lon, Real &lat);
void XYPFromRLL(Real lon, Real lat, Real &dX, Real &dY, int &nP);

}  // namespace CubedSphereUtility

#endif  // SRC_EXO3_CUBED_SPHERE_UTILITY_HPP_
