#ifndef SRC_EXO3_WRAPED_COORDINATE_HPP_
#define SRC_EXO3_WRAPED_COORDINATE_HPP_

// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>

class MeshBlock;
class ParameterInput;


class WarpedCoordinate : public Coordinates {
 public:
  WarpedCoordinate(MeshBlock *pmb, ParameterInput *pin, bool flag);

  void Face1Area(const int k, const int j, const int il, const int iu,
                 AthenaArray<Real> &area) final;
  void Face2Area(const int k, const int j, const int il, const int iu,
                 AthenaArray<Real> &area) final;
  void Face3Area(const int k, const int j, const int il, const int iu,
                 AthenaArray<Real> &area) final;
  Real GetFace1Area(const int k, const int j, const int i) final;
  Real GetFace2Area(const int k, const int j, const int i) final;
  Real GetFace3Area(const int k, const int j, const int i) final;

  void CellVolume(const int k, const int j, const int il, const int iu,
                  AthenaArray<Real> &vol);
  Real GetCellVolume(const int k, const int j, const int i);

  void AddCoordTermsDivergence(const Real dt, const AthenaArray<Real> *flux,
                               const AthenaArray<Real> &prim,
                               const AthenaArray<Real> &bcc,
                               AthenaArray<Real> &u) final;
  private:
    Real w_x1, c0, c1, c2;
    Real w2(Real x1);
    Real int_w2(Real x1);
    Real int_w2x1(Real x1);
    Real dw2(Real x1l, Real x1r);
    Real mean_w2(Real x1l, Real x1r);
};

#endif  // SRC_EXO3_WRAPED_COORDINATE_HPP_
