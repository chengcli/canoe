#ifndef SRC_EXO3_WRAPED_COORDINATE_HPP_
#define SRC_EXO3_WRAPED_COORDINATE_HPP_

// athena
#include <athena/athena.hpp>
#include <athena/coordinates/coordinates.hpp>

class MeshBlock;
class ParameterInput;

# define w2_c0 1
# define w2_c1 1e-3

template <typename R> R w2(R x1) {
    return w2_c0 + w2_c1 * x1;
}

template <typename R> R dw2(R x1l, R x1r) {
    return w2_c1;
}

template <typename R> R mean_w2(R x1l, R x1r) {
    return w2_c0 + w2_c1 * ((x1l + x1r) / 2);
}

template <typename R> R mean_w2x1(R x1l, R x1r) {
    return (w2_c0 * ((x1l + x1r) / 2)
      + w2_c1 * ((x1l*x1l + x1l*x1r + x1r*x1r) / 3));
}

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
};

#endif  // SRC_EXO3_WRAPED_COORDINATE_HPP_
