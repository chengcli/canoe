#ifndef SRC_EXO3_GNOMONIC_EQUIANGLE_HPP_
#define SRC_EXO3_GNOMONIC_EQUIANGLE_HPP_

// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>

class MeshBlock;
class ParameterInput;

class GnomonicEquiangle : public Coordinates {
 public:
  GnomonicEquiangle(MeshBlock *pmb, ParameterInput *pin, bool flag);
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
                  AthenaArray<Real> &vol) final;
  Real GetCellVolume(const int k, const int j, const int i) final;

  void CenterWidth1(const int k, const int j, const int il, const int iu,
                    AthenaArray<Real> &dx1) final;
  void CenterWidth2(const int k, const int j, const int il, const int iu,
                    AthenaArray<Real> &dx2) final;
  void CenterWidth3(const int k, const int j, const int il, const int iu,
                    AthenaArray<Real> &dx3) final;

  void CellMetric(const int k, const int j, const int il, const int iu,
                  AthenaArray<Real> &g, AthenaArray<Real> &g_inv) final;
  void Face1Metric(const int k, const int j, const int il, const int iu,
                   AthenaArray<Real> &g, AthenaArray<Real> &g_inv) final;
  void Face2Metric(const int k, const int j, const int il, const int iu,
                   AthenaArray<Real> &g, AthenaArray<Real> &g_inv) final;
  void Face3Metric(const int k, const int j, const int il, const int iu,
                   AthenaArray<Real> &g, AthenaArray<Real> &g_inv) final;

  void PrimToLocal1(const int k, const int j, const int il, const int iu,
                    const AthenaArray<Real> &b1_vals,
                    AthenaArray<Real> &prim_left, AthenaArray<Real> &prim_right,
                    AthenaArray<Real> &bx) final;
  void PrimToLocal2(const int k, const int j, const int il, const int iu,
                    const AthenaArray<Real> &b2_vals,
                    AthenaArray<Real> &prim_left, AthenaArray<Real> &prim_right,
                    AthenaArray<Real> &bx) final;
  void PrimToLocal3(const int k, const int j, const int il, const int iu,
                    const AthenaArray<Real> &b3_vals,
                    AthenaArray<Real> &prim_left, AthenaArray<Real> &prim_right,
                    AthenaArray<Real> &bx) final;

  void FluxToGlobal1(const int k, const int j, const int il, const int iu,
                     const AthenaArray<Real> &cons,
                     const AthenaArray<Real> &bbx, AthenaArray<Real> &flux,
                     AthenaArray<Real> &ey, AthenaArray<Real> &ez) final;
  void FluxToGlobal2(const int k, const int j, const int il, const int iu,
                     const AthenaArray<Real> &cons,
                     const AthenaArray<Real> &bbx, AthenaArray<Real> &flux,
                     AthenaArray<Real> &ey, AthenaArray<Real> &ez) final;
  void FluxToGlobal3(const int k, const int j, const int il, const int iu,
                     const AthenaArray<Real> &cons,
                     const AthenaArray<Real> &bbx, AthenaArray<Real> &flux,
                     AthenaArray<Real> &ey, AthenaArray<Real> &ez) final;
  void AddCoordTermsDivergence(const Real dt, const AthenaArray<Real> *flux,
                               const AthenaArray<Real> &prim,
                               const AthenaArray<Real> &bcc,
                               AthenaArray<Real> &u) final;

  Real GetCosineCell(const int k, const int j) const {
    return cosine_cell_kj_(k, j);
  }

  Real GetSineCell(const int k, const int j) const {
    return sine_cell_kj_(k, j);
  }

 protected:
  Real sphericalTri(Real x1, Real x2, Real x3, Real y1, Real y2, Real y3);

  AthenaArray<Real> cosine_cell_kj_, sine_cell_kj_;
  AthenaArray<Real> cosine_face2_kj_, sine_face2_kj_;
  AthenaArray<Real> cosine_face3_kj_, sine_face3_kj_;
  AthenaArray<Real> x_ov_rD_kji_, y_ov_rC_kji_;
  AthenaArray<Real> dx2f_ang_kj_, dx3f_ang_kj_;
  AthenaArray<Real> dx2f_ang_face3_kj_, dx3f_ang_face2_kj_;
};

#endif  // SRC_EXO3_GNOMONIC_EQUIANGLE_HPP_
