#ifndef SRC_EXO3_AFFINE_COORDINATE_HPP_
#define SRC_EXO3_AFFINE_COORDINATE_HPP_

// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>

class MeshBlock;
class ParameterInput;

class AffineCoordinate : public Coordinates {
 public:
  AffineCoordinate(MeshBlock *pmb, ParameterInput *pin, bool flag);

  void Face1Area(const int k, const int j, const int il, const int iu,
                 AthenaArray<Real> &area) final;
  void Face2Area(const int k, const int j, const int il, const int iu,
                 AthenaArray<Real> &area) final;
  void Face3Area(const int k, const int j, const int il, const int iu,
                 AthenaArray<Real> &area) final;
  Real GetFace1Area(const int k, const int j, const int i) final;
  Real GetFace2Area(const int k, const int j, const int i) final;
  Real GetFace3Area(const int k, const int j, const int i) final;

  void VolCenterFace1Area(const int k, const int j, const int il, const int iu,
                          AthenaArray<Real> &area) final;
  void VolCenterFace2Area(const int k, const int j, const int il, const int iu,
                          AthenaArray<Real> &area) final;
  void VolCenterFace3Area(const int k, const int j, const int il, const int iu,
                          AthenaArray<Real> &area) final;
  void CellVolume(const int k, const int j, const int il, const int iu,
                  AthenaArray<Real> &vol);
  Real GetCellVolume(const int k, const int j, const int i);

  void CellMetric(const int k, const int j, const int il, const int iu,
                  AthenaArray<Real> &g, AthenaArray<Real> &g_inv);
  void Face1Metric(const int k, const int j, const int il, const int iu,
                   AthenaArray<Real> &g, AthenaArray<Real> &g_inv);
  void Face2Metric(const int k, const int j, const int il, const int iu,
                   AthenaArray<Real> &g, AthenaArray<Real> &g_inv);
  void Face3Metric(const int k, const int j, const int il, const int iu,
                   AthenaArray<Real> &g, AthenaArray<Real> &g_inv);

  void PrimToLocal2(const int k, const int j, const int il, const int iu,
                    const AthenaArray<Real> &b1_vals,
                    AthenaArray<Real> &prim_left, AthenaArray<Real> &prim_right,
                    AthenaArray<Real> &bx);
  void PrimToLocal3(const int k, const int j, const int il, const int iu,
                    const AthenaArray<Real> &b1_vals,
                    AthenaArray<Real> &prim_left, AthenaArray<Real> &prim_right,
                    AthenaArray<Real> &bx);

  void FluxToGlobal2(const int k, const int j, const int il, const int iu,
                     const AthenaArray<Real> &cons,
                     const AthenaArray<Real> &bbx, AthenaArray<Real> &flux,
                     AthenaArray<Real> &ey, AthenaArray<Real> &ez);
  void FluxToGlobal3(const int k, const int j, const int il, const int iu,
                     const AthenaArray<Real> &cons,
                     const AthenaArray<Real> &bbx, AthenaArray<Real> &flux,
                     AthenaArray<Real> &ey, AthenaArray<Real> &ez);

 protected:
  Real theta_;
  Real sin_theta_;
  Real cos_theta_;
};

#endif  // SRC_EXO3_AFFINE_COORDINATE_HPP_
